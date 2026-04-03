#!/bin/bash
set -euo pipefail

###############################################################################
# Unified r/m estimation pipeline
#
# This script runs the full workflow:
#  1) Prodigal gene prediction
#  2) OrthoFinder single-copy ortholog detection
#  3) Codon-based core-gene supermatrix construction (MAFFT + pal2nal + AMAS)
#  4) IQ-TREE phylogeny + unrooting with ape (R)
#  5) ClonalFrameML to estimate R/theta, delta, nu, and r/m
#
# USAGE:
#   1) Edit the USER SETTINGS section below.
#   2) Make executable:  chmod +x rm_full_pipeline.sh
#   3) Run:             ./rm_full_pipeline.sh
###############################################################################

######################################
# ==== USER SETTINGS (EDIT HERE) ====
######################################

# Directory containing your input MAG assemblies (FASTA format)
ASSEMBLY_DIR="/path/to/My_Project/1_input_assemblies"

# Root directory where all results will be written
OUTPUT_ROOT="/path/to/My_Project/3_analysis_results"

# Conda environment with: prodigal, orthofinder, mafft, pal2nal, AMAS.py, iqtree/iqtree2, R+ape, python+biopython
BIOINFO_ENV_NAME="r_m_pipeline"

# Conda environment with: clonalframeml
CFML_ENV_NAME="clonalframeml_env"

# Number of threads to use for MAFFT / OrthoFinder / IQ-TREE
N_THREADS=32

######################################
# ==== NO CHANGES NORMALLY BELOW ====
######################################

echo "==============================="
echo "  r/m full pipeline started"
echo "==============================="
echo "ASSEMBLY_DIR : ${ASSEMBLY_DIR}"
echo "OUTPUT_ROOT  : ${OUTPUT_ROOT}"
echo "BIOINFO_ENV  : ${BIOINFO_ENV_NAME}"
echo "CFML_ENV     : ${CFML_ENV_NAME}"
echo "THREADS      : ${N_THREADS}"
echo

if [ ! -d "${ASSEMBLY_DIR}" ]; then
  echo "ERROR: ASSEMBLY_DIR does not exist: ${ASSEMBLY_DIR}" >&2
  exit 1
fi

mkdir -p "${OUTPUT_ROOT}"

# Load conda
if command -v conda >/dev/null 2>&1; then
  # shellcheck disable=SC1091
  source "$(conda info --base)/etc/profile.d/conda.sh"
else
  echo "ERROR: conda command not found. Please load Conda before running this script." >&2
  exit 1
fi

input_folder_name="$(basename "${ASSEMBLY_DIR}")"

PRODIGAL_DIR="${OUTPUT_ROOT}/${input_folder_name}_prodigal_o"
ORTHOF_DIR="${OUTPUT_ROOT}/${input_folder_name}_orthofinder_o"
SUPERMATRIX_DIR="${OUTPUT_ROOT}/${input_folder_name}_supermatrix_o"
SCRIPTS_DIR="${OUTPUT_ROOT}/${input_folder_name}_scripts"
PY_SCRIPT="${SCRIPTS_DIR}/create_codon_supermatrix.py"

mkdir -p "${PRODIGAL_DIR}" "${ORTHOF_DIR}" "${SUPERMATRIX_DIR}" "${SCRIPTS_DIR}"

echo "Working directories:"
echo "  PRODIGAL    : ${PRODIGAL_DIR}"
echo "  ORTHOFINDER : ${ORTHOF_DIR}"
echo "  SUPERMATRIX : ${SUPERMATRIX_DIR}"
echo "  SCRIPTS     : ${SCRIPTS_DIR}"
echo

##########################################
# Phase 0: write helper Python script
##########################################

echo "--- Phase 0: Writing helper Python script (create_codon_supermatrix.py) ---"

cat << 'PYEOF' > "${PY_SCRIPT}"
import os
import sys
import subprocess
import argparse
from collections import defaultdict

from Bio import SeqIO


def is_alignment_valid(filepath: str) -> bool:
    # Return True if FASTA file has at least one non-empty sequence.
    try:
        for record in SeqIO.parse(filepath, "fasta"):
            if len(record.seq) > 0:
                return True
    except Exception:
        return False
    return False


def main(orthofinder_results_dir: str, cds_dir: str, output_dir: str, threads: int) -> None:
    single_copy_ogs_file = os.path.join(
        orthofinder_results_dir, "Orthogroups", "Orthogroups_SingleCopyOrthologues.txt"
    )
    if not os.path.isfile(single_copy_ogs_file):
        print(f"ERROR: Cannot find {single_copy_ogs_file}", file=sys.stderr)
        sys.exit(1)

    with open(single_copy_ogs_file) as f:
        single_copy_ogs = {line.strip() for line in f if line.strip()}

    print(f"{len(single_copy_ogs)} single-copy orthologous groups found.")

    orthogroups_tsv_file = os.path.join(
        orthofinder_results_dir, "Orthogroups", "Orthogroups.tsv"
    )
    if not os.path.isfile(orthogroups_tsv_file):
        print(f"ERROR: Cannot find {orthogroups_tsv_file}", file=sys.stderr)
        sys.exit(1)

    og_to_genes = defaultdict(list)
    with open(orthogroups_tsv_file) as f:
        header = f.readline().rstrip("\n").split("\t")
        species_list = [s.replace(".faa", "") for s in header[1:]]

        for line in f:
            parts = line.rstrip("\n").split("\t")
            og = parts[0]
            if og not in single_copy_ogs:
                continue
            for i, gene_list_str in enumerate(parts[1:]):
                genes = gene_list_str.split(", ")
                if genes and genes[0]:
                    og_to_genes[og].append((species_list[i], genes[0]))

    codon_alignments_dir = os.path.join(output_dir, "codon_alignments")
    temp_files_dir = os.path.join(output_dir, "temp_files")
    os.makedirs(codon_alignments_dir, exist_ok=True)
    os.makedirs(temp_files_dir, exist_ok=True)

    valid_aligned_files = []

    print("Creating per-orthogroup codon alignments...")
    for og, gene_info_list in og_to_genes.items():
        unaligned_cds_file = os.path.join(temp_files_dir, f"{og}_cds.fna")
        unaligned_pep_file = os.path.join(temp_files_dir, f"{og}_pep.faa")
        aligned_pep_file = os.path.join(temp_files_dir, f"{og}_pep.aln")
        final_codon_aln_file = os.path.join(codon_alignments_dir, f"{og}_codon.aln")

        with open(unaligned_cds_file, "w") as f_cds, open(unaligned_pep_file, "w") as f_pep:
            for species_name, long_gene_id in gene_info_list:
                clean_gene_id = long_gene_id.split(" ")[0]
                cds_file_path = os.path.join(cds_dir, f"{species_name}.fna")
                if not os.path.isfile(cds_file_path):
                    print(f"WARNING: CDS file not found for {species_name}: {cds_file_path}", file=sys.stderr)
                    continue
                found = False
                for record in SeqIO.parse(cds_file_path, "fasta"):
                    if record.id == clean_gene_id:
                        # write CDS
                        record_for_cds = record[:]  # shallow copy
                        record_for_cds.id = species_name
                        record_for_cds.description = ""
                        SeqIO.write(record_for_cds, f_cds, "fasta")
                        # write AA (translated)
                        aa_seq = record.seq.translate(table=11, to_stop=True)
                        record_for_aa = record[:]
                        record_for_aa.seq = aa_seq
                        record_for_aa.id = species_name
                        record_for_aa.description = ""
                        SeqIO.write(record_for_aa, f_pep, "fasta")
                        found = True
                        break
                if not found:
                    print(
                        f"WARNING: gene {clean_gene_id} not found in {cds_file_path}",
                        file=sys.stderr,
                    )

        # Align protein sequences with MAFFT
        mafft_cmd = [
            "mafft",
            "--auto",
            "--thread",
            str(threads),
            unaligned_pep_file,
        ]
        with open(aligned_pep_file, "w") as out_f:
            subprocess.run(
                mafft_cmd,
                check=True,
                stdout=out_f,
                stderr=subprocess.DEVNULL,
            )

        # Convert to codon alignment using pal2nal
        pal2nal_cmd = [
            "pal2nal.pl",
            aligned_pep_file,
            unaligned_cds_file,
            "-output",
            "fasta",
            "-nogap",
        ]
        with open(final_codon_aln_file, "w") as out_f:
            subprocess.run(
                pal2nal_cmd,
                check=True,
                stdout=out_f,
                stderr=subprocess.DEVNULL,
            )

        if is_alignment_valid(final_codon_aln_file):
            valid_aligned_files.append(final_codon_aln_file)
        else:
            print(
                f"WARNING: {og} codon alignment has no valid sequences; excluded from concatenation.",
                file=sys.stderr,
            )

    print(f"\n{len(valid_aligned_files)} valid codon alignments will be concatenated.")
    if not valid_aligned_files:
        print("ERROR: No valid codon alignments found. Aborting.", file=sys.stderr)
        sys.exit(1)

    output_supermatrix_path = os.path.join(output_dir, "my_codon_supermatrix.fasta")
    amas_cmd = [
        "AMAS.py",
        "concat",
        "-f",
        "fasta",
        "-d",
        "dna",
        "-i",
    ] + valid_aligned_files + [
        "-t",
        output_supermatrix_path,
    ]
    subprocess.run(amas_cmd, check=True)

    print(f"\nDone! Concatenated codon alignment written to: {output_supermatrix_path}")

    # clean up temp directory
    import shutil

    shutil.rmtree(temp_files_dir, ignore_errors=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Create codon-based supermatrix from OrthoFinder results.")
    parser.add_argument("--orthofinder_dir", required=True, help="Path to OrthoFinder Results_* directory")
    parser.add_argument("--cds_dir", required=True, help="Path to directory with Prodigal CDS .fna files")
    parser.add_argument("--output_dir", required=True, help="Output directory for supermatrix and intermediate files")
    parser.add_argument("--threads", type=int, default=8, help="Threads for MAFFT")
    args = parser.parse_args()
    main(
        orthofinder_results_dir=args.orthofinder_dir,
        cds_dir=args.cds_dir,
        output_dir=args.output_dir,
        threads=args.threads,
    )
PYEOF

echo "  -> Python helper written to: ${PY_SCRIPT}"
echo

##########################################
# Phase 1: Prodigal gene prediction
##########################################

echo "--- Phase 1: Prodigal gene prediction ---"

FAA_OUTPUT_DIR="${PRODIGAL_DIR}/faa_files"
FNA_OUTPUT_DIR="${PRODIGAL_DIR}/fna_files"
mkdir -p "${FAA_OUTPUT_DIR}" "${FNA_OUTPUT_DIR}"

echo "Input assemblies : ${ASSEMBLY_DIR}"
echo "Output (FAA)     : ${FAA_OUTPUT_DIR}"
echo "Output (FNA)     : ${FNA_OUTPUT_DIR}"
echo

conda activate "${BIOINFO_ENV_NAME}"

shopt -s nullglob
assembly_files=( "${ASSEMBLY_DIR}"/*.fna "${ASSEMBLY_DIR}"/*.fa "${ASSEMBLY_DIR}"/*.fasta "${ASSEMBLY_DIR}"/*.contigs )
if [ ${#assembly_files[@]} -eq 0 ]; then
  echo "ERROR: No assembly files (*.fna|*.fa|*.fasta|*.contigs) found in ${ASSEMBLY_DIR}" >&2
  exit 1
fi

for assembly_file in "${assembly_files[@]}"; do
  base_filename="$(basename "${assembly_file}")"
  basename_noext="${base_filename%.*}"
  echo "  - Running Prodigal on: ${base_filename}"
  prodigal -i "${assembly_file}" \
           -a "${FAA_OUTPUT_DIR}/${basename_noext}.faa" \
           -d "${FNA_OUTPUT_DIR}/${basename_noext}.fna"
done

echo "--- Phase 1 completed ---"
echo

##########################################
# Phase 2: OrthoFinder
##########################################

echo "--- Phase 2: OrthoFinder ---"
echo "Input (FAA)  : ${FAA_OUTPUT_DIR}"
echo "Output root  : ${ORTHOF_DIR}"
echo

orthofinder -f "${FAA_OUTPUT_DIR}" -o "${ORTHOF_DIR}" -t "${N_THREADS}"

echo "--- Phase 2 completed ---"
echo

##########################################
# Phase 3: Codon supermatrix
##########################################

echo "--- Phase 3: Codon-based supermatrix construction ---"

ORTHOF_RESULTS_SUBDIR="$(find "${ORTHOF_DIR}" -type d -name 'Results_*' | head -n 1 || true)"
if [ -z "${ORTHOF_RESULTS_SUBDIR}" ]; then
  echo "ERROR: Could not find Orthofinder Results_* directory under ${ORTHOF_DIR}" >&2
  exit 1
fi


echo "OrthoFinder results: ${ORTHOF_RESULTS_SUBDIR}"
echo "CDS directory      : ${FNA_OUTPUT_DIR}"
echo "Supermatrix out    : ${SUPERMATRIX_DIR}"
echo

python "${PY_SCRIPT}" \
  --orthofinder_dir "${ORTHOF_RESULTS_SUBDIR}" \
  --cds_dir "${FNA_OUTPUT_DIR}" \
  --output_dir "${SUPERMATRIX_DIR}" \
  --threads "${N_THREADS}"


echo "--- Phase 3 completed ---"
echo

##########################################
# Phase 4: IQ-TREE + ape + ClonalFrameML
##########################################

echo "--- Phase 4: IQ-TREE phylogeny & ClonalFrameML ---"

cd "${SUPERMATRIX_DIR}"
echo "Working directory: $(pwd)"
echo

SUPERMATRIX_FASTA="my_codon_supermatrix.fasta"
TREEFILE="${SUPERMATRIX_FASTA}.treefile"
UNROOTED_TREE="final.unrooted.tre"
CFML_PREFIX="my_final_clonalframe_analysis"

if [ ! -f "${SUPERMATRIX_FASTA}" ]; then
  echo "ERROR: Supermatrix FASTA not found: ${SUPERMATRIX_FASTA}" >&2
  exit 1
fi


echo "[4-1] IQ-TREE maximum-likelihood tree + 1000 ultrafast bootstrap"
if [ -f "${TREEFILE}" ]; then
  echo "  - Treefile already exists (${TREEFILE}). Skipping IQ-TREE."
else
  iqtree -s "${SUPERMATRIX_FASTA}" -B 1000 -T AUTO
fi


echo "[4-2] Unrooting tree with R/ape -> ${UNROOTED_TREE}"
Rscript -e "library(ape); tree <- read.tree('${TREEFILE}'); unrooted <- unroot(tree); write.tree(unrooted, file='${UNROOTED_TREE}');"


echo "[4-3] Switching to ClonalFrameML environment: ${CFML_ENV_NAME}"
conda activate "${CFML_ENV_NAME}"


echo "  - Running ClonalFrameML..."
ClonalFrameML "${UNROOTED_TREE}" "${SUPERMATRIX_FASTA}" "${CFML_PREFIX}"


echo

echo "======================================="
echo "  ALL ANALYSIS STAGES COMPLETED"
echo "  Final ClonalFrameML prefix:"
echo "    $(pwd)/${CFML_PREFIX}"
echo "======================================="

