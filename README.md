# Individual-Project
Guide to reproducing report results

This document provides step-by-step instructions to reproduce the quantitative result reported in this work: the protein backbone and ligand pose RMSD (relative to the experimental structure) for the CMP-sialic acid transporter bound to 5-methyl CMP (PDB 6XBO), including the reproducibility statistic (mean &plusmn; SD across three independent pipeline runs).

<img width="2818" height="1663" alt="image" src="https://github.com/user-attachments/assets/77a00f0f-4420-4a2f-961e-c0ce781c42ec" />

## Software available note
This work used a pre-release build of the CCD2MD package which includes features not yet available via the standard `pip install CCD2MD`installation method documented by the package. Because this build cannot be currently installed independently, a full snapshot of the exact package is provided as `ccd2md_testing_snapshot.tar.gz` in this repository, capturing both code and reference data.

To install, create a conda/mamba environment with GROMACS 2025.4:
>conda create -n my_env python=3.12 -y\
conda activate my_env\
conda install -c conda-forge gromacs=2025.4\
pip install GromacsWrapper numpy pandas scipy h5py matplotlib jax requests\
pip install "MemPrOD==0.0.9"\
pip install "cg2at-lite==0.3.1"\
pip install vermouth

>tar -xzf ccd2md_testing_snapshot.tar.gz\
cp -r ccd2md/ /path/to/new/env/lib/python3.12/site-packages

## Step 1: Generate the AlphaFold3 prediction
1. Build an environment fold AlphaFold3 using the installation guide: https://github.com/google-deepmind/alphafold3/blob/main/docs/installation.md
2. Usw the provided .json file as input (6XBO.json):
  Format:
```
  {
  "dialect": "alphafold3",
  "version": 2,
  "name": "6XBO-AF3",
  "sequences": [
    {
      "protein": {
        "id": "A",
        "sequence": "MAPARENVSLFFKLYCLTVMTLVAAAYTVALRYTRTTAEELYFSTTAVCITEVIKLLISVGLLAKETGSLGRFKASLSENVLGSPKELAKLSVPSLVYAVQNNMAFLALSNLDAAVYQVTYQLKIPCTALCTVLMLNRTLSKLQWISVFMLCGGVTLVQWKPAQATKVVVAQNPLLGFGAIAIAVLCSGFAGVYFEKVLKSSDTSLWVRNIQMYLSGIVVTLAGTYLSDGAEIQEKGFFYGYTYYVWFVIFLASVGGLYTSVVVKYTDNIMKGFSAAAAIVLSTIASVLLFGLQITLSFALGALLVCVSIYLYGLPRQDTTSNSLEVLFQ"
      }
    },
    {
      "ligand": {
        "id": "B",
        "ccdCodes": [
          "5MC"
        ]
      }
    }
  ],
  "modelSeeds": [
    1
  ],
  "bondedAtomPairs": null,
  "userCCDPath": "/storage/chem/lfsmgr/SRG/Packages/CCD2MD/CCD2MD.cif"
}
```

&rarr; This includes the sequence of the protein and the CCD code for the ligand (5MC)

Run:
```
python run_alphafold.py --json_path 6XBO.json --output_dir <output_dir>
```
3. Output: `6XBO-AF3_model.cif` (to be used as input for the automated pipeline in the next step)

## Step 2: Run the automated pipeline

```
chmod +x automate_script.sh
```

```
./automate_script.sh -i 6XBO-AF3_model.cif -n 6XBO-5MC -g 5MC -Mem
```

Pipeline outline:
<img width="988" height="325" alt="image" src="https://github.com/user-attachments/assets/b735636e-b30c-4210-9bd1-304bfb3cfeb5" />

This produces `final_system.pdb` and `topol.top` in the working directory, representing the membrane-embedded, ligand-reinserted atomistic system.

## Step 3: Energy minimisation of the final system
A standard steepest-descent energy minimisation `em.mdp` is sufficient. The specific `em.mdp` used for the reported result is provided alongside this guide. 

```
gmx grompp -f em.mdp -c final_system.pdb -p topol.top -o em.tpr
```
```
gmx mdrun -deffnm em -v -ntmpi 1
```
This produces `em.gro`

## Step 4: Compute RMSD 
First prepare the crystal structure and energy minimised structure for `gmx rms` which compares atom-for-atom between a reference structure (-s) and a trajectory (-f). It requires both to have the same number of atoms in the same order. A crystal PDB and an MD/AF3-derived .gro may not satisfy that directly due to offset residue numbering, different atom names, the .gro having solvent/ions/lipids/hydrogens while the crystal does not.
```
python prepare_for_gmx_rms.py --crystal 6XBO.pdb --ligand 5MC --pred em.gro
```
```
echo -e "0\n0" | gmx rms -s crystal_ref.pdb -f pred_1.pdb -n index.ndx -o protein_rmsd.xvg
```
```
echo -e "0\n2" | gmx rms -s crystal_ref.pdb -f pred_1.pdb -n index.ndx -o ligand.xvg
```
Both `.xvg` files report the RMSD, in nm, in the second column of their single data row.

## Expected result
To reproduce the 'reproducible result', repeat steps 2-4 twice more from the same input file (`6XBO-AF3_model.cif`), using a fresh working directory for each run. 

Collect the three protein and three ligand RMSD values, then compute the mean and standard deviation.

Compute mean:
```
python -c "
import statistics as stats
protein = [<run1>, <run2>, <run3>]   # values in Angstrom
ligand  = [<run1>, <run2>, <run3>]
print(f'Protein: {stats.mean(protein):.3f} +/- {stats.stdev(protein):.3f} A')
print(f'Ligand:  {stats.mean(ligand):.3f} +/- {stats.stdev(ligand):.3f} A')
"
```
The reported results are

`Protein: 4.086 $\pm$ 0.010 Å`

`Ligand: 1.371 $\pm$ 0.202 Å`
