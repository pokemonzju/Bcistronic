#  `extract_transcript_orf_phylocsf_score.py`

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



# `calculate_reads_coverage_base_bamfile.py`

## Overview

This script analyzes the alignment of reads from a BAM file against annotated ORFs (Open Reading Frames) provided in a CSV file. It calculates how well each read covers the ORFs and summarizes statistics across all reads, including average coverage scores, match ratios, and coverage patterns. We use **avg_ORF_match_score** as the **ORF-normalized_read_coverage_score** in article.


## Input Files

- **BAM file**: Contains aligned sequencing reads (e.g., RNA-seq, Ribo-seq) in binary alignment format.
- **ORF CSV file**: Describes genomic locations of ORFs using the following structure:
    - `Chromosome`: Chromosome name
    - `Strand`: '+' or '-'
    - `ORF_genome_location`: e.g., `ORF1_[(100-200); (300-400)] | ORF2_[(500-600)]`


## Features

1. **Parse ORF file**: Extract ORF names and segment intervals.
2. **Build interval trees**: Create genomic interval maps for fast overlap checking.
3. **Parse BAM file**:
   - For each read, determine overlaps with ORF regions.
   - Compute per-read ORF match score.
   - Track which ORFs were covered and the extent of coverage.
4. **Summarize results**:
   - Average per-read match score
   - Average match ratio for each ORF
   - Number of reads that match only a single ORF
   - Number of reads that match multiple ORFs
   - Formatted summary strings for ORF match information

## Output Files

- `[trans_id]_result.csv`: Detailed results per read, including:
    - `Reads_name`
    - `Reads_ORF_score`
    - `Covered_ORFs`
    - `ORF_match_info`
    - `ORF_match_ratio`

- `[trans_id]_merge.csv`: Summary statistics including:
    - `avg_Reads_ORF_score`
    - `avg_ORF_match_score`
    - `avg_ORF_match_info`
    - `avg_ORF_match_ratio`
    - `Single_ORF_reads`
    - `Positive_ORF_reads`


## Example Usage

```bash
python calculate_reads_coverage_base_bamfile.py \
  --bam_file sample.bam \
  --orf_file transcript123.csv
```

## Requirements

- Python package
- `pysam`
- `pandas`
- `intervaltree`


