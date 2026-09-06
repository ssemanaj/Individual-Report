#!/bin/bash
set -euo pipefail

# script for preparing system post-AF3
echo "starting"
input_cif=""
base_name=""
ligand=""
maxwarn=2
membrane=false
upper_leaflet=()
lower_leaflet=()
posi_ion="NA"
nega_ion="CL"
conc=0.15

usage() {
    cat <<EOF
Usage: $0 -i <input.cif> -n <base_name> -g <ligand_resname> [-w maxwarn] [-Mem] [-u LIPID:RATIO ...] [-l LIPID:RATIO ...] [-C conc] [-pname posi_ion] [-nname nega_ion]

    -i     Input AF3 .cif file (required)
    -n     Base name for intermediate/output files, e.g. "KcsA-TEA" (required)
    -g     Ligand residue name to exclude from membrane building / reinsert (required)
    -w     maxwarn value for ccd2at (default: $maxwarn)
    -Mem   Enable membrane embedding (-mem flag for ccd2at)
    -u     Upper leaflet lipid:ratio, e.g. -u POPE:7, -u POPG:2 (default: pure POPC)
    -l     Lower leaflet lipid:ratio, e.g. -l POPE:7, -l POPG:2 (default: pure POPC)
    -C     Ion concentration (default: $conc)
    -pname Positive ion name (default: $posi_ion)
    -nname Negative ion name (default: $nega_ion)
    -h     Show this help
EOF
    exit 1
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        -i) input_cif="$2"; shift 2 ;;
        -n) base_name="$2"; shift 2 ;;
        -g) ligand="$2"; shift 2 ;;
        -w) maxwarn="$2"; shift 2 ;;
        -Mem) membrane=true; shift ;;
        -u) upper_leaflet+=("$2"); shift 2 ;;
        -l) lower_leaflet+=("$2"); shift 2 ;;
        -C) conc="$2"; shift 2 ;;
        -pname) posi_ion="$2"; shift 2 ;;
        -nname) nega_ion="$2"; shift 2 ;;
        -h) usage ;;
        *) echo "Uknown option: $1" >&2; usage ;;
    esac
done

#while getopts "i:n:g:w:Memu:l:C:pname:nname:h" option; do
#    case $option in
#        i) input_cif="$OPTARG" ;;
#        n) base_name="$OPTARG" ;;
#        g) ligand="$OPTARG" ;;
#        w) maxwarn="$OPTARG" ;;
#        Mem) membrane=true ;;
#        u) upper_leaflet+=("$OPTARG") ;;
#        l) lower_leaflet+=("$OPTARG") ;;
#        C) conc="$OPTARG" ;;
#        pname) posi_ion="$OPTARG" ;;
#        nname) nega_ion="$OPTARG" ;;
#        h) usage ;;
#        \?) usage ;;
#        :) echo "Option -$OPTARG requires an argument" >&2; usage ;;
#    esac
#done

[ -z "$input_cif" ] && { echo "Error: -i input.cif is required" >&2; usage; }
[ -z "$base_name" ] && { echo "Error: -n base_name is required" >&2; usage; }
[ -z "$ligand" ] && { echo "Error: -g ligand_name is required" >&2; usage; }
[ -f "$input_cif" ] || { echo "Error: input file '$input_cif' not found" >&2; exit 1; }

detect_ions() {
    local cif_file="$1"
    local known_ions=("NA" "K" "CL" "MG" "CA" "ZN" "FE2" "FE" "MN" "CO" "NI" "CU1" "CU")

    # Pull unique values from the _atom_site.label_comp_id (or auth_comp_id) column only
    local comp_ids
    comp_ids=$(awk '
        /^loop_/ { col=0; compcol=0; next }
        /^_atom_site\./ {
            col++
            if ($0 ~ /_atom_site\.label_comp_id/ || $0 ~ /_atom_site\.auth_comp_id/) {
                compcol = col
            }
            next
        }
        compcol > 0 && NF >= compcol && ($1 == "ATOM" || $1 == "HETATM") {
            print $compcol
        }
    ' "$cif_file" | sort -u)

    detected_ions=()
    for ion in "${known_ions[@]}"; do
        if echo "$comp_ids" | grep -qxw "$ion"; then
            detected_ions+=("$ion")
        fi
    done
}

# 1. ccd2at for system preparation (membrane insertion if membrane protein, waterbox for non-membrane)
aligned_ligand_pdb="${base_name}-ligand-fitted.pdb"

mem_args=()
ccd2at_args=()

if [ "$membrane" = true ]; then
# membrane case
    echo -e "STEP 1: Building membrane system (excluding ligand: $ligand)\n"

    detect_ions "$input_cif"

    ignore_targets=("$ligand")
    if [ ${#detected_ions[@]} -gt 0 ]; then
        echo "# INFO: Detected ions in structure: ${detected_ions[*]}"
        ignore_targets+=("${detected_ions[@]}")
    else
        echo "# INFO: No ions detected."
    fi

    ccd2at_membrane="${base_name}_membrane.pdb"
    ccd2at_ligand_only="${base_name}_only.pdb"
    ccd2at_ligand_only_H="${base_name}_only_H.pdb"

    #mem_args+=("-mem")
    for lipid in "${upper_leaflet[@]}"; do
        mem_args+=("-u" "$lipid")
    done
    for lipid in "${lower_leaflet[@]}"; do
        mem_args+=("-l" "$lipid")
    done
    mem_args+=("-posi_c0" "${posi_ion}:1")
    mem_args+=("-negi_c0" "${nega_ion}:1")
    
    echo -e "STEP 1A: processing system with no membrane\n"
    mkdir ./protein_ligand_only
    cp "$input_cif" ./protein_ligand_only
    cd ./protein_ligand_only
    ccd2at_hdf5 "$input_cif" "$ccd2at_ligand_only" -St
    cd ../
    echo -e "STEP 1A complete: output $ccd2at_ligand_only_H and topology in directory ./protein_ligand_only/\n"

    echo -e "STEP 1B: processing protein only with membrane\n"
    mkdir ./protein_membrane
    cp "$input_cif" ./protein_membrane
    cd ./protein_membrane
    ccd2at_hdf5 "$input_cif" "$ccd2at_membrane" -C "$conc" -il -M -maxwarn 1 -mem "${mem_args[@]}"
    cd ../
    echo -e "STEP 1B complete: output $ccd2at_membrane and topology in ./protein_membrane/\n"

    echo -e "STEP 2: inserting ligands/ions into membrane system.\n"
    echo -e "4\n4" | gmx confrms -f1 "./protein_membrane/${ccd2at_membrane}" -f2 "./protein_ligand_only/${ccd2at_ligand_only_H}" -one -o "$aligned_ligand_pdb"      # backbone for aligning groups
    grep -v -E "END" "./protein_membrane/$ccd2at_membrane" > temporary_system.pdb
    for target in "${ignore_targets[@]}"; do
        grep "$target" "$aligned_ligand_pdb" | grep -v "TITLE" >> temporary_system.pdb
    done
    echo "ENDMDL" >> temporary_system.pdb
    mv temporary_system.pdb final_system.pdb

    mkdir ./topol/
    cp ./protein_membrane/*.itp ./topol/
    cp ./protein_membrane/topol.top ./
    cp ./protein_ligand_only/posre_Other*.itp ./topol/
    cp ./protein_ligand_only/topol_Other*.itp ./topol/
    cp -r protein_ligand_only/charmm36-ccd2md.ff ./topol/

    sed -i "/^\[ system \]/i #include \"topol_Other.itp\"" "topol.top"

    if grep -q "^#ifdef DISRES" "topol.top"; then
        sed -i "/^#ifdef DISRES/i Other    1" "topol.top"
    else
        echo "Other    1" >> "topol.top"
    fi

    echo "# INFO: Merged 'Other' (count 1) from topol_Other.itp"

    sed -i 's|#include "|#include "./topol/|' topol.top
    
else
#non-membrane case
    echo "Step 1: Building waterbox system"
    ccd2at_output="${base_name}_box.pdb"
    ccd2at_args+=("-I" "-pname" "${posi_ion}" "-nname" "${nega_ion}" "-neutral")
    ccd2at_hdf5 "$input_cif" "$ccd2at_output" -C "$conc" "${ccd2at_args[@]}"
    
    mkdir ./topol/
    mv *.itp topol/
    sed -i 's|#include "|#include "./topol/|' topol.top
fi

echo "Done."
