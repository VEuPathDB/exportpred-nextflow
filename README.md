# exportpred-nextflow

Nextflow pipeline that predicts protein export signals with ExportPred and publishes
the results as a genome browser track.

## Overview

This pipeline is part of VEuPathDB's genomic data workflows. It runs the ExportPred
tool over a proteome to identify N-terminal export/host-targeting signals — the
short motifs (e.g. PEXEL-type `RLE`/`KLD` classes) that direct proteins for export
into the host cell, most notably in apicomplexan parasites such as *Plasmodium*. For
each predicted signal, ExportPred's raw output is parsed into a GFF3 feature that
decomposes the N-terminal region into its annotated segments (initial methionine,
signal/leader peptide, hydrophobic region, and export motif), scored by the
underlying model. The merged results are sorted, bgzip-compressed, and tabix-indexed
for direct use as a genome/protein browser track.

## Requirements

- [Nextflow](https://www.nextflow.io/) (DSL2)
- A container engine: Docker (default) or Singularity/Apptainer

Containers used:
- `veupathdb/exportpred:1.0.0` — the ExportPred prediction tool
- `bioperl/bioperl:stable` — conversion of ExportPred output to GFF3
- `biocontainers/tabix:v1.9-11-deb_cv1` — sorting, bgzip, and tabix indexing

## Usage

```
nextflow run VEuPathDB/exportpred-nextflow -r main -resume \
  --inputFilePath /path/to/proteins.fa \
  --outputDir /path/to/output \
  -C conf/docker.config
```

To run under Singularity/Apptainer (e.g. on an HPC cluster via LSF):

```
nextflow run VEuPathDB/exportpred-nextflow -r main -resume \
  --inputFilePath /path/to/proteins.fa \
  --outputDir /path/to/output \
  -C conf/lsf.config
```

There is a single, unnamed entry point — no `-entry` flag is needed. The input
FASTA is split into subsets (`fastaSubsetSize`) that are scored by ExportPred in
parallel (up to `process.maxForks = 5` concurrent tasks).

## Key Parameters

| Parameter | Default | Description |
|---|---|---|
| `inputFilePath` | `data/input.fa` | FASTA file of protein sequences to score |
| `fastaSubsetSize` | `500` | Number of sequences per subset chunk processed in each parallel ExportPred task |
| `outputFileName` | `exportPred.out` | Base name of the merged GFF output before compression |
| `outputDir` | `$launchDir/output` | Directory where output files are published |

## Output

Written to `outputDir`:

- **`exportPred.out.gz`** (+ `.gz.tbi`) — bgzip-compressed, tabix-indexed GFF3 of
  predicted export signals across the input proteome. Each protein with a predicted
  signal has a parent feature (type `RLE` or `KLD`, scored) spanning the motif, plus
  child sub-features marking the decomposed regions of the signal (e.g. `a-met`,
  `a-leader`, `a-hydrophobic`).
