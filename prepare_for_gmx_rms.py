"""
Prepare matched protein + ligand structures from a crystal PDB and one or more
predicted GROMACS .gro structures, so that native `gmx rms` can be used to
compute RMSD (instead of a custom Kabsch implementation).

writes out:
  - crystal_ref.pdb   : the crystal atoms we could confidently match, in a
                       fixed, canonical atom order
  - pred_<n>.pdb  : the same atoms (same order, same count) from each
                       predicted .gro file
  - index.ndx        : GROMACS index groups "Protein_CA", "Protein_Backbone",
                       and "Ligand", pointing at that shared atom order

After that, plain `gmx rms` (see the printed commands at the end) does the
actual least-squares fitting + RMSD calculation -- no custom alignment code
needed for that part.

HOW ATOMS ARE MATCHED
----------------------
Protein: each chain's residue sequence (one-letter code) is aligned between
the PDB and the GRO chain with difflib.SequenceMatcher (the GRO chain is
assumed to be the same sequence, just potentially longer / gap-free). Only
residues present in BOTH are kept.

Ligand: PDB and GRO atom names for the ligand often differ, so atoms are
matched by their bonded graph (heavy atoms only, bonds = pairs within 1.75 A)
using networkx graph isomorphism on element type. Where the ligand has a
genuine topological symmetry (e.g. equivalent phosphate oxygens), the
isomorphism that minimizes the raw (unfitted) RMSD is kept, so the ordering
is picked once, consistently, before writing any files.

USAGE
-----
    python3 prepare_for_gmx_rms.py --crystal CRYSTAL.pdb --ligand LIGAND_RESNAME --pred PREDICTED.gro [PREDICTED2.gro ...]
    python3 prepare_for_gmx_rms.py --crystal CRYSTAL.pdb --ligand LIGAND_RESNAME --pred PREDICTED.gro --chains A B --outdir results/

Chains default to every chain ID that appears in the crystal PDB's ATOM
records (in the order they first appear), so for a single-chain structure
you don't need to pass --chains at all.
"""
import argparse
import difflib
import os
import numpy as np
import networkx as nx

PROT_RES = set("ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL".split())
THREE_TO_ONE = {
    'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C','GLN':'Q','GLU':'E','GLY':'G',
    'HIS':'H','ILE':'I','LEU':'L','LYS':'K','MET':'M','PHE':'F','PRO':'P','SER':'S',
    'THR':'T','TRP':'W','TYR':'Y','VAL':'V'}
BB_ATOMS = ["N", "CA", "C", "O"]


def detect_chains(pdb_path):
    """Return protein chain IDs in the order they first appear in ATOM records."""
    order = []
    seen = set()
    with open(pdb_path) as f:
        for line in f:
            if line.startswith("ATOM"):
                c = line[21]
                if c not in seen:
                    seen.add(c)
                    order.append(c)
    return order


def parse_pdb(path, chains, ligand_resname):
    prot_seq = {c: [] for c in chains}          # chain -> [(resid, resname), ...]
    prot_atoms = {}                              # (chain, resid, atomname) -> coord
    lig_atoms = {}                                # atomname -> coord (assumes ONE ligand copy)
    prev = None
    with open(path) as f:
        for line in f:
            if not (line.startswith("ATOM") or line.startswith("HETATM")):
                continue
            if line[16] not in (' ', 'A'):        # skip alt-locs other than blank/A
                continue
            chain = line[21]
            resid = int(line[22:26])
            resname = line[17:20].strip()
            atomname = line[12:16].strip()
            coord = np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])])
            if line.startswith("ATOM") and chain in chains:
                key = (chain, resid)
                if key != prev:
                    prot_seq[chain].append((resid, resname))
                    prev = key
                prot_atoms[(chain, resid, atomname)] = coord
            elif line.startswith("HETATM") and resname == ligand_resname:
                lig_atoms[atomname] = coord
    return prot_seq, prot_atoms, lig_atoms


def parse_gro(path, ligand_resname):
    with open(path) as f:
        lines = f.readlines()
    natoms = int(lines[1])
    body = lines[2:2 + natoms]
    segments = []      # list of [(resid, resname), ...] per chain, in file order
    prot_atoms = {}     # (chain_idx, resid, atomname) -> coord
    lig_atoms = {}       # atomname -> coord
    chain_idx, cur = -1, []
    prev_resid, prev_was_prot = None, False
    for l in body:
        resname = l[5:10].strip()
        resid = int(l[0:5])
        atomname = l[10:15].strip()
        coord = np.array([float(l[20:28]), float(l[28:36]), float(l[36:44])]) * 10.0  # nm -> A
        is_prot = resname in PROT_RES
        if is_prot:
            if not prev_was_prot or (prev_resid is not None and resid == 0 and prev_resid != 0):
                if cur:
                    segments.append(cur)
                cur = [(resid, resname)]
                chain_idx += 1
            elif resid != cur[-1][0]:
                cur.append((resid, resname))
            prot_atoms[(chain_idx, resid, atomname)] = coord
        else:
            if prev_was_prot and cur:
                segments.append(cur)
                cur = []
            if resname == ligand_resname:
                lig_atoms[atomname] = coord
        prev_resid, prev_was_prot = resid, is_prot
    if cur:
        segments.append(cur)
    return segments, prot_atoms, lig_atoms


def align_chain(pdb_list, gro_list):
    """Return {pdb_resid: gro_resid} for residues common to both, via sequence alignment."""
    pdb_str = ''.join(THREE_TO_ONE.get(r, 'X') for _, r in pdb_list)
    gro_str = ''.join(THREE_TO_ONE.get(r, 'X') for _, r in gro_list)
    sm = difflib.SequenceMatcher(None, gro_str, pdb_str, autojunk=False)
    mapping = {}
    for block in sm.get_matching_blocks():
        for i in range(block.size):
            pdb_resid = pdb_list[block.b + i][0]
            gro_resid = gro_list[block.a + i][0]
            mapping[pdb_resid] = gro_resid
    return mapping


def build_heavy_graph(atomdict, cutoff=1.75):
    names = [n for n in atomdict if not n.startswith('H')]
    coords = np.array([atomdict[n] for n in names])
    G = nx.Graph()
    for n in names:
        G.add_node(n, elem=n[0])
    for i in range(len(names)):
        for j in range(i + 1, len(names)):
            if np.linalg.norm(coords[i] - coords[j]) < cutoff:
                G.add_edge(names[i], names[j])
    return G


def match_ligand_atoms(pdb_lig, gro_lig):
    """Return list of (pdb_atomname, gro_atomname) pairs, resolving any
    topological symmetry by minimum raw (unfitted) distance."""
    GP, GG = build_heavy_graph(pdb_lig), build_heavy_graph(gro_lig)
    GM = nx.algorithms.isomorphism.GraphMatcher(GP, GG, node_match=lambda a, b: a['elem'] == b['elem'])
    best, best_cost = None, None
    for mapping in GM.isomorphisms_iter():
        cost = sum(np.linalg.norm(pdb_lig[k] - gro_lig[v]) for k, v in mapping.items())
        if best_cost is None or cost < best_cost:
            best_cost, best = cost, mapping
    if best is None:
        raise RuntimeError("Ligand heavy-atom graphs are not isomorphic -- check ligand identity/naming.")
    return list(best.items())


def write_pdb(path, atoms, chain_id="A"):
    """atoms: list of (resname, resid, atomname, coord, record_type).
    Writes strict PDB fixed-column format so downstream tools (GROMACS, etc.)
    parse atom name / resname / chain / resid correctly:
      cols 1-6 record, 7-11 serial, 13-16 atom name, 17 altLoc (blank),
      18-20 resName, 22 chainID, 23-26 resSeq, 31-54 x/y/z, 55-66 occ/B.
    """
    with open(path, "w") as f:
        for i, (resname, resid, atomname, coord, rec) in enumerate(atoms, start=1):
            # PDB convention: atom names <4 chars get a leading blank (cols 13-16)
            name_field = f"{atomname:<4}" if len(atomname) >= 4 else f" {atomname:<3}"
            f.write(
                f"{rec:<6}{i:>5} {name_field}{' '}{resname:>3} {chain_id}{resid:>4}    "
                f"{coord[0]:8.3f}{coord[1]:8.3f}{coord[2]:8.3f}{1.0:6.2f}{0.0:6.2f}\n"
            )
        f.write("END\n")


def parse_args():
    p = argparse.ArgumentParser(
        description="Prepare matched protein+ligand structures from a crystal PDB and "
                     "one or more predicted .gro files, for use with `gmx rms`.")
    p.add_argument("--crystal", required=True, help="Crystal structure, e.g. 6XBO.pdb")
    p.add_argument("--ligand", required=True, dest="ligand_resname",
                    help="Three/four-letter ligand residue name, e.g. 5MC, SLB, PHO")
    p.add_argument("--pred", required=True, nargs="+", dest="gro",
                    help="One or more predicted .gro files (replicas of the same system)")
    p.add_argument("--chains", nargs="+", default=None,
                    help="PDB chain IDs to include, in order (default: auto-detect all chains in the PDB)")
    p.add_argument("--outdir", default="./", help="Directory to write output files to")
    return p.parse_args()


def main():
    args = parse_args()
    pdb_path = args.crystal
    gro_paths = args.gro
    ligand_resname = args.ligand_resname
    chains = args.chains if args.chains else detect_chains(pdb_path)
    outdir = args.outdir

    os.makedirs(outdir, exist_ok=True)
    print(f"Using chains: {chains}")

    pdb_prot_seq, pdb_prot_atoms, pdb_lig = parse_pdb(pdb_path, chains, ligand_resname)

    # Use the FIRST gro file just to build the shared chain/residue mapping and
    # canonical atom ordering; later gro files are assumed to share the same
    # topology/atom order (true for replicate MD runs of the same system).
    gro_segments0, _, _ = parse_gro(gro_paths[0], ligand_resname)
    if len(gro_segments0) != len(chains):
        raise RuntimeError(f"Expected {len(chains)} protein chain(s) in {gro_paths[0]}, found {len(gro_segments0)}")

    resid_maps = {}   # chain -> {pdb_resid: gro_resid}
    for c, seg in zip(chains, gro_segments0):
        resid_maps[c] = align_chain(pdb_prot_seq[c], seg)

    # gro's own protein atoms (from the first replica) are needed here too, to
    # skip backbone atoms that aren't present under the same name in gro (e.g.
    # C-terminal residues often use OT1/OT2 instead of O for the carbonyl O).
    _, gro_prot_atoms0, _ = parse_gro(gro_paths[0], ligand_resname)

    # Canonical ordered list of (chain, pdb_resid, atomname) for backbone atoms
    backbone_order = []
    ca_flags = []       # parallel list: True if this entry is a CA atom
    for c in chains:
        chain_idx = chains.index(c)
        for pdb_resid, _ in pdb_prot_seq[c]:
            if pdb_resid not in resid_maps[c]:
                continue
            gro_resid = resid_maps[c][pdb_resid]
            for a in BB_ATOMS:
                if (c, pdb_resid, a) in pdb_prot_atoms and (chain_idx, gro_resid, a) in gro_prot_atoms0:
                    backbone_order.append((c, pdb_resid, a))
                    ca_flags.append(a == "CA")

    # Ligand atom order, resolved once against the first gro file
    _, _, gro_lig0 = parse_gro(gro_paths[0], ligand_resname)
    ligand_pairs = match_ligand_atoms(pdb_lig, gro_lig0)   # [(pdb_name, gro_name), ...]

    # actual amino-acid resname per (chain, resid), so the written PDB carries
    # real residue names (helps GROMACS's own default "Protein" group work too,
    # as a fallback if -n index.ndx is forgotten)
    resname_lookup = {(c, r): rn for c in chains for r, rn in pdb_prot_seq[c]}

    # ---- write reference (crystal) structure ----
    ref_atoms = []
    for (c, resid, aname) in backbone_order:
        ref_atoms.append((resname_lookup[(c, resid)], resid, aname, pdb_prot_atoms[(c, resid, aname)], "ATOM"))
    for pdb_name, _ in ligand_pairs:
        ref_atoms.append((ligand_resname, 9999, pdb_name, pdb_lig[pdb_name], "HETATM"))
    write_pdb(f"{outdir}/crystal_ref.pdb", ref_atoms)

    # ---- write each mobile (predicted) structure, same atom order ----
    mobile_paths = []
    for gi, gpath in enumerate(gro_paths, start=1):
        segments, prot_atoms, lig_atoms = parse_gro(gpath, ligand_resname)
        mob_atoms = []
        for idx, (c, resid, aname) in enumerate(backbone_order):
            chain_idx = chains.index(c)
            gro_resid = resid_maps[c][resid]
            mob_atoms.append((resname_lookup[(c, resid)], resid, aname, prot_atoms[(chain_idx, gro_resid, aname)], "ATOM"))
        for pdb_name, gro_name in ligand_pairs:
            mob_atoms.append((ligand_resname, 9999, pdb_name, lig_atoms[gro_name], "HETATM"))
        outpath = f"{outdir}/pred_{gi}.pdb"
        write_pdb(outpath, mob_atoms)
        mobile_paths.append(outpath)

    # ---- write GROMACS index file ----
    n_bb = len(backbone_order)
    n_lig = len(ligand_pairs)
    ca_idx = [i + 1 for i, flag in enumerate(ca_flags) if flag]          # 1-based
    bb_idx = list(range(1, n_bb + 1))
    lig_idx = list(range(n_bb + 1, n_bb + n_lig + 1))

    def write_ndx_group(f, name, idx_list):
        f.write(f"[ {name} ]\n")
        for i in range(0, len(idx_list), 15):
            f.write(" ".join(str(x) for x in idx_list[i:i + 15]) + "\n")

    with open(f"{outdir}/index.ndx", "w") as f:
        write_ndx_group(f, "Protein_CA", ca_idx)
        write_ndx_group(f, "Protein_Backbone", bb_idx)
        write_ndx_group(f, "Ligand", lig_idx)

    print(f"Wrote {outdir}/crystal_ref.pdb  ({n_bb} protein backbone atoms incl. {len(ca_idx)} CA, "
          f"{n_lig} ligand heavy atoms)")
    for p in mobile_paths:
        print(f"Wrote {p}")
    print(f"Wrote {outdir}/index.ndx  (groups: Protein_CA, Protein_Backbone, Ligand)")

    print("\n--- Next steps: run these with GROMACS ---\n")
    if len(mobile_paths) > 1:
        print("# 1. Concatenate replicas into one pseudo-trajectory (one 'frame' per replica)")
        print(f"gmx trjcat -f {' '.join(mobile_paths)} -o {outdir}/mobile_traj.xtc -cat -settime <<EOF")
        for i in range(len(mobile_paths)):
            print("0")   # accept default start time offsets non-interactively
        print("EOF")
        traj = f"{outdir}/mobile_traj.xtc"
    else:
        traj = mobile_paths[0]

    print(f"\n# 2. Protein RMSD (fit on CA, calculate on CA) after fit rot+trans")
    print(f"gmx rms -s {outdir}/crystal_ref.pdb -f {traj} -n {outdir}/index.ndx "
          f"-o {outdir}/protein_ca_rmsd.xvg -fit rot+trans <<EOF")
    print("Protein_CA")
    print("Protein_CA")
    print("EOF")

    print(f"\n# 3. Protein backbone RMSD (fit + calc on N,CA,C,O)")
    print(f"gmx rms -s {outdir}/crystal_ref.pdb -f {traj} -n {outdir}/index.ndx "
          f"-o {outdir}/protein_bb_rmsd.xvg -fit rot+trans <<EOF")
    print("Protein_Backbone")
    print("Protein_Backbone")
    print("EOF")

    print(f"\n# 4. Ligand POSE RMSD: fit on protein CA, but calculate RMSD on the ligand")
    print(f"#    (this is the standard 'docking pose' RMSD -- how far the ligand sits")
    print(f"#     from its crystal position once the protein itself is aligned)")
    print(f"gmx rms -s {outdir}/crystal_ref.pdb -f {traj} -n {outdir}/index.ndx "
          f"-o {outdir}/ligand_pose_rmsd.xvg -fit rot+trans <<EOF")
    print("Protein_CA")
    print("Ligand")
    print("EOF")

    print(f"\n# 5. Mean +/- SD across replicas (if you concatenated multiple), from any .xvg:")
    print(f"python3 -c \"import numpy as np; d=np.loadtxt('{outdir}/ligand_pose_rmsd.xvg', comments=('#','@')); "
          f"print(d[:,1].mean(), d[:,1].std(ddof=1))\"")


if __name__ == "__main__":
    main()
