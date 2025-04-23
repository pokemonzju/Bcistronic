#  `extract_transcript_orf_phylocsff_score.py`

## Overview

This script is designed to **map ORF transcript coordinates to genomic coordinates** and **evaluate coding potential using PhyloCSF scores**. 

## Features

The script performs the following steps:

1. **Parse the GFF3 file** to extract exon structures of transcripts.
2. **Load the ORF CSV file**, which includes start and end coordinates of ORFs.
3. **Map ORF transcript positions to genomic coordinates**.
4. Extract ORF start locations and calculate:
   - Relative start differences
   - Modulo 3 remainder of start positions
5. Adjust:
   - Genomic start positions for correct frame identification
   - Translation frame annotations (e.g., frame+1, frame-2)
6. **Extract and analyze ORF end positions** similarly.
7. **Calculate corrected genomic intervals**.
8. **Determine WIG file paths** for each ORF segment.
9. Convert coordinates to **bedGraph format**.
10. **Query BigWig files for PhyloCSF scores** using UCSC's `bigWigSummary`.
11. Compute **average PhyloCSF scores per ORF** and assess whether all scores are positive.

## Input Files

| Argument | Description |
|----------|-------------|
| `--gff3_file` | GFF3 annotation file with transcript and exon information |
| `--orf_csv_file` | CSV file containing `Trans_ID` and `ORF_ID`, where `ORF_ID` follows format like `ORF3_ENST00000335137.7:100:400` |
| `--base_directory` | Base directory containing per-chromosome PhyloCSF BigWig files like `output_chr1/PhyloCSFRaw+1.bw` |
| `--output_file` | Path to the final output CSV file |

## Output

The output CSV includes:

- `Chromosome`: Chromosome containing the ORF
- `Strand`: Strand orientation
- `adj_ORF_genome_location`: Adjusted genomic coordinates of ORFs
- `frame_info`: Reading frame of each ORF segment
- `phylocsf_file`: BigWig file paths used for score retrieval
- `ORF_phylocsf_score`: PhyloCSF score for each ORF segment
- `ORF_phylocsf_merge_score`: Average PhyloCSF score per ORF
- `ORF_all_positive`: Whether all ORF segments have positive PhyloCSF scores

## Example Usage

python extract_transcript_orf_phylocsff_score.py \
  --gff3_file transcripts.gff3 \
  --orf_csv_file orfs.csv \
  --output_file results.csv \
  --base_directory /data/phylocsf_tracks/

## Requirements

- Python 3.6+
- `pandas`
- UCSC's `bigWigSummary` tool

