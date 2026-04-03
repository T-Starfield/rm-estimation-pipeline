# rm-estimation-pipeline

Pipeline for estimating recombination-to-mutation ratios (r/m) in microbial genomes using single-copy core genes.

This workflow was developed for the analysis of metagenome-assembled genomes (MAGs) from subseafloor sediments.

---

## Overview

The pipeline performs:

1. Protein sequence alignment (MAFFT)  
2. Codon-aware nucleotide alignment (pal2nal)  
3. Concatenation of gene alignments  
4. Phylogenetic tree-based analysis  
5. Estimation of recombination parameters using ClonalFrameML  

Outputs include R/theta, 1/delta, nu, and derived r/m values.

---

## Requirements

The following software must be installed and available in your PATH:

- MAFFT (v7 or later)
- pal2nal (v14 or later)
- ClonalFrameML

Optional:
- IQ-TREE (for phylogenetic tree inference)

---

## Input

The pipeline expects:

- Amino acid sequences of single-copy core genes (FASTA format)
- Corresponding nucleotide sequences (FASTA format)
- Consistent genome/sample identifiers across files

---

## Usage

```bash
bash rm_full_pipeline.sh <input_directory> <output_directory>

---

## Output

The pipeline generates:
- Codon alignments for each gen
- Concatenated alignment file
- ClonalFrameML output:
	-- R/theta
	-- 1/delta
	-- nu
	-- r/m estimates
		

