#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pysam
import csv
import argparse
import re
import os
import pandas as pd
from intervaltree import Interval, IntervalTree
from collections import defaultdict, OrderedDict

# ---------------------- Step 1: Step 1: Parse ORF data ----------------------

def parse_orf_file(orf_file):
    """Parse the ORF file"""
    orf_data = defaultdict(lambda: defaultdict(lambda: {'orfs': [], 'transcript_orf_count': 0}))
    
    with open(orf_file, 'r', newline='', encoding='utf-8') as csvfile:
        reader = csv.DictReader(csvfile, delimiter=',')
        for row in reader:
            chromosome = row['Chromosome']
            strand = row['Strand']
            orf_genome_location = row['ORF_genome_location']
            orf_genome_entries = orf_genome_location.split('|')

            orfs = []
            for orf_genome_entry in orf_genome_entries:
                parts = orf_genome_entry.strip().split('_[')
                if len(parts) != 2:
                    continue
                orf_name = parts[0]
                segments = [tuple(map(int, re.findall(r'\d+', seg))) for seg in parts[1].rstrip(']').split(';')]
                orfs.append({'name': orf_name, 'segments': segments})

            orf_data[chromosome][strand]['orfs'].extend(orfs)
            orf_data[chromosome][strand]['transcript_orf_count'] = len(orfs)

    return orf_data

def build_interval_trees(orf_data):
    """Build ORF IntervalTree"""
    orf_trees = defaultdict(lambda: defaultdict(IntervalTree))
    for chrom in orf_data:
        for strand in orf_data[chrom]:
            for orf in orf_data[chrom][strand]['orfs']:
                for segment in orf['segments']:
                    orf_trees[chrom][strand][segment[0]:segment[1]] = orf['name']
    return orf_trees

# ---------------------- Step 2: Step 2: Parse BAM and compute ORF matching ----------------------

def parse_bam_and_compute_scores(bam_file, orf_data, orf_trees, trans_id):
    """Step 2: Parse BAM and compute ORF matching率"""
    samfile = pysam.AlignmentFile(bam_file, "rb")
    reads_with_scores = []
    chrom, strand = list(orf_data.keys())[0], list(orf_data[list(orf_data.keys())[0]].keys())[0]
    orf_lengths = {orf['name']: sum(e - s for s, e in orf['segments']) for orf in orf_data[chrom][strand]['orfs']}

    for read in samfile:
        read_name = read.query_name
        regions = [(read.reference_start, read.reference_start + read.reference_length)]
        orf_overlap_info = defaultdict(int)

        for s, e in regions:
            for iv in orf_trees[chrom][strand][s:e]:
                orf_name = iv.data
                overlap_length = min(e, iv.end) - max(s, iv.begin)
                orf_overlap_info[orf_name] += overlap_length

        reads_orf_score = sum(orf_overlap_info.values()) / sum(orf_lengths.values()) if orf_overlap_info else 0.0
        orf_match_info = " | ".join(f"{orf}_[{overlap}/{orf_lengths[orf]}]" for orf, overlap in orf_overlap_info.items())
        orf_match_ratio = " | ".join(f"{orf}_[{overlap/orf_lengths[orf]:.2f}]" for orf, overlap in orf_overlap_info.items())

        reads_with_scores.append([trans_id, read_name, reads_orf_score, "; ".join(orf_overlap_info.keys()), orf_match_info, orf_match_ratio])

    samfile.close()
    return reads_with_scores

# ---------------------- Step 3: Step 3: Compute ORF statistics ----------------------

def extract_avg_orf_match_info(orf_match_infos):
    """Compute average match value from ORF_match_info"""
    orf_match_counts = defaultdict(list)
    orf_lengths = {}
    
    for match_info in orf_match_infos.dropna():
        for entry in match_info.split(" | "):
            match = re.search(r'(ORF\d+)\_\[(\d+)/(\d+)\]', entry)
            if match:
                orf_name, match_value, orf_length = match.groups()
                orf_match_counts[orf_name].append(int(match_value))
                orf_lengths[orf_name] = orf_length

    return " | ".join(f"{orf}_[{sum(matches) / len(matches):.2f}/{orf_lengths[orf]}]" 
                      for orf, matches in orf_match_counts.items())

def process_orf_results(reads_with_scores, trans_id, transcript_orf_count):
    """Process ORF match results and generate summary statistics"""
    df = pd.DataFrame(reads_with_scores, columns=["Trans_ID", "Reads_name", "Reads_ORF_score", "Covered_ORFs", "ORF_match_info", "ORF_match_ratio"])

    avg_orf_match_scores = df["ORF_match_ratio"].apply(
        lambda x: sum(map(float, re.findall(r'\d+\.\d+', str(x)))) / transcript_orf_count
        if pd.notna(x) and transcript_orf_count > 0 else 0
    ).mean()

    def extract_avg_ratio_per_orf(orf_match_ratios):
        """Compute average match ratio per ORF"""
        orf_scores = defaultdict(list)
        orf_order = []

        for ratios in orf_match_ratios.dropna():
            if isinstance(ratios, str):
                for entry in ratios.split(" | "):
                    match = re.search(r'(ORF\d+)\_\[(\d+\.\d+)\]', entry)
                    if match:
                        orf_name, score = match.groups()
                        orf_scores[orf_name].append(float(score))
                        if orf_name not in orf_order:
                            orf_order.append(orf_name)

        avg_ratios = OrderedDict()
        for orf in orf_order:
            if orf in orf_scores:
                avg_ratios[orf] = sum(orf_scores[orf]) / len(orf_scores[orf])

        return " | ".join(f"{orf}_[{avg_ratios[orf]:.2f}]" for orf in avg_ratios.keys())

    def count_single_orf_reads(covered_orfs):
        """Count reads that cover only a single ORF"""
        orf_counts = defaultdict(int)

        for orfs in covered_orfs.dropna():
            orf_list = [o.strip() for o in orfs.split("; ") if o.strip()]
            if len(orf_list) == 1:
                orf_counts[orf_list[0]] += 1

        return " | ".join(f"{orf}_[{count}]" for orf, count in orf_counts.items() if orf)

    def count_positive_orf_reads(covered_orfs, total_reads):
        """Count reads that cover multiple ORFs"""
        multi_orf_reads = sum(1 for orfs in covered_orfs.dropna() if len([o.strip() for o in orfs.split("; ") if o.strip()]) >= 2)
        return f"[{multi_orf_reads}/{total_reads}]"

    output_df = pd.DataFrame({
        "Trans_ID": [trans_id],
        "avg_Reads_ORF_score": [df["Reads_ORF_score"].mean()],
        "avg_ORF_match_score": [avg_orf_match_scores],
        "avg_ORF_match_info": [extract_avg_orf_match_info(df["ORF_match_info"])],
        "avg_ORF_match_ratio": [extract_avg_ratio_per_orf(df["ORF_match_ratio"])],
        "Single_ORF_reads": [count_single_orf_reads(df["Covered_ORFs"])],
        "Positive_ORF_reads": [count_positive_orf_reads(df["Covered_ORFs"], len(df))]
    })

    return df, output_df

# ---------------------- Step 4: Step 4: Main function ----------------------

def main():
    parser = argparse.ArgumentParser(description="Parse BAM and ORF data and compute ORF match statistics")
    parser.add_argument("--bam_file", required=True, help="Input BAM file")
    parser.add_argument("--orf_file", required=True, help="Input ORF CSV file")
    args = parser.parse_args()

    orf_data = parse_orf_file(args.orf_file)
    orf_trees = build_interval_trees(orf_data)
    trans_id = os.path.splitext(os.path.basename(args.orf_file))[0]

    reads_with_scores = parse_bam_and_compute_scores(args.bam_file, orf_data, orf_trees, trans_id)
    df_result, df_merge = process_orf_results(reads_with_scores, trans_id, orf_data[list(orf_data.keys())[0]][list(orf_data[list(orf_data.keys())[0]].keys())[0]]['transcript_orf_count'])

    df_result.to_csv(f"{trans_id}_result.csv", index=False)
    df_merge.to_csv(f"{trans_id}_merge.csv", index=False)

if __name__ == "__main__":
    main()
