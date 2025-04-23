#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# usage: extract_transcript_orf_phylocsff_score_correct_frame.py [-h] --gff3_file GFF3_FILE --orf_csv_file ORF_CSV_FILE --output_file OUTPUT_FILE --base_directory BASE_DIRECTORY

import csv
import re
from collections import defaultdict
import pandas as pd
import sys
import os
import subprocess
from statistics import mean

# -----------------------------
# Function definitions for Steps 1 to 7
# -----------------------------

def parse_gff(gff_file):
    """
    Parse the GFF3 file to extract exon information for each transcript.
    """
    transcripts = defaultdict(list)
    
    with open(gff_file, 'r') as f:
        for line in f:
            if line.startswith("#"):
                continue
            
            parts = line.strip().split('\t')
            if len(parts) < 9:
                continue
            
            chrom, source, feature_type, start, end, score, strand, phase, attributes = parts
            start = int(start)
            end = int(end)
            
            # Extract transcript_id
            transcript_id = None
            for attr in attributes.split(";"):
                if attr.startswith("transcript_id="):
                    transcript_id = attr.split("=")[1]
                    break
            
            # Collect only exon information
            if transcript_id and feature_type == "exon":
                transcripts[transcript_id].append((chrom, start, end, strand))
    
    # Sort exons of each transcript by start position
    for exons in transcripts.values():
        exons.sort(key=lambda x: x[1])
    
    return transcripts

def parse_orf_csv(orf_csv_file):
    """
    Parse the ORF CSV file, extract ORF data, and retain the original columns.
    """
    print(f"Parsing ORF CSV file: {orf_csv_file}")
    orf_data = {}
    original_data = []
    with open(orf_csv_file) as f:
        reader = csv.DictReader(f)
        for row in reader:
            trans_id = row['Trans_ID']
            orf_info = row['ORF_ID'].split(';')
            orfs = []
            for orf in orf_info:
                try:
                    # Split ORF information, e.g., ORF5_ENST00000003100.13:2482:2691
                    if "_" in orf and ":" in orf:
                        orf_name, rest = orf.split('_', 1)
                        transcript_id, start, end = rest.split(':')
                        start, end = int(start)+1, int(end)+1
                        orfs.append({
                            "ORF Name": orf_name,
                            "Transcript ID": transcript_id,
                            "Start": start,
                            "End": end
                        })
                    else:
                        print(f"Warning: ORF '{orf}' format is unexpected, skipping.")
                except ValueError as e:
                    print(f"Error parsing ORF '{orf}': {e}")
            # Store multiple ORF entries for each transcript
            orf_data[trans_id] = orfs
            # Retain the original data
            original_data.append(row)
    print(f"Parsed {len(orf_data)} transcripts with ORF data.")
    return orf_data, original_data

def build_position_map(exons, strand):
    """
    Create a mapping from transcript coordinates to genomic coordinates for each transcript.
    """
    position_map = []
    transcript_pos = 1
    if strand == "+":
        for exon_start, exon_end in exons:
            for genome_pos in range(exon_start, exon_end + 1):
                position_map.append((transcript_pos, genome_pos))
                transcript_pos += 1
    else:
        for exon_start, exon_end in reversed(exons):
            for genome_pos in range(exon_end, exon_start - 1, -1):
                position_map.append((transcript_pos, genome_pos))
                transcript_pos += 1
    return position_map

def map_orf_to_genome(transcript_exons, orf_data):
    """
    Map ORF transcript coordinates to genomic coordinates.
    """
    results = []

    for trans_id, orfs in orf_data.items():
        if trans_id not in transcript_exons:
            print(f"Warning: Transcript {trans_id} not found in GFF3 data")
            continue
        
        exons = [(start, end) for _, start, end, _ in transcript_exons[trans_id]]
        chrom = transcript_exons[trans_id][0][0]  # Extract chromosome number
        strand = transcript_exons[trans_id][0][3]  # Extract strand information
        
        # Construct the position_map of the current transcript
        position_map = build_position_map(exons, strand)
        
        # Map each ORF individually
        orf_locations = []
        orf_trans_locations = []
        for orf in orfs:
            orf_id = orf["ORF Name"]
            orf_start = orf["Start"]
            orf_end = orf["End"]
            genome_coords = []
            trans_coords = []
            
            for trans_pos, genome_pos in position_map:
                if orf_start <= trans_pos <= orf_end:
                    genome_coords.append(genome_pos)
                    trans_coords.append(trans_pos)
            
            # Merge consecutive coordinates into intervals
            if genome_coords:
                genome_ranges = []
                trans_ranges = []
                
                range_start = genome_coords[0]
                range_end = genome_coords[0]
                trans_range_start = trans_coords[0]
                trans_range_end = trans_coords[0]
                
                for i in range(1, len(genome_coords)):
                    # Determine if they are consecutive
                    if (strand == "+" and genome_coords[i] == range_end + 1) or (strand == "-" and genome_coords[i] == range_end - 1):
                        range_end = genome_coords[i]
                        trans_range_end = trans_coords[i]
                    else:
                        # Save the previous interval
                        if strand == "+":
                            genome_ranges.append((range_start, range_end))
                            trans_ranges.append((trans_range_start, trans_range_end))
                        else:
                            # The negative strand needs to be reversed
                            genome_ranges.append((range_end, range_start))
                            trans_ranges.append((trans_range_end, trans_range_start))
                        
                        range_start = genome_coords[i]
                        range_end = genome_coords[i]
                        trans_range_start = trans_coords[i]
                        trans_range_end = trans_coords[i]
                
                # Process the last interval
                if strand == "+":
                    genome_ranges.append((range_start, range_end))
                    trans_ranges.append((trans_range_start, trans_range_end))
                else:
                    genome_ranges.append((range_end, range_start))
                    trans_ranges.append((trans_range_end, trans_range_start))
                
                # Format
                orf_location_str = f"{orf_id}_[{'; '.join(f'({s}-{e})' for s, e in genome_ranges)}]"
                orf_trans_location_str = f"{orf_id}_[{'; '.join(f'({s}-{e})' for s, e in trans_ranges)}]"
                orf_locations.append(orf_location_str)
                orf_trans_locations.append(orf_trans_location_str)
        
        # merge all ORF
        results.append({
            "Trans_ID": trans_id,
            "Chromosome": chrom,
            "Strand": strand,
            "ORF_genome_location": " | ".join(orf_locations),
            "ORF_trans_location": " | ".join(orf_trans_locations)
        })
    
    return results

def extract_orf_start_locations(results):
    """
    Extract the start position information of ORFs from ORF_trans_location and ORF_genome_location.
    """
    for result in results:
        trans_location = result["ORF_trans_location"]
        orf_trans_start_locations = []
        for orf in trans_location.split(" | "):
            trans_matches = re.findall(r"\((\d+)-(\d+)\)", orf)
            if trans_matches:
                start_positions = [f"({min(int(m[0]), int(m[1]))})" for m in trans_matches]
                orf_trans_start_locations.append(f"{orf.split('_')[0]}_[{'; '.join(start_positions)}]")
            else:
                orf_trans_start_locations.append(f"{orf.split('_')[0]}_[(NA)]")
        result["ORF_trans_start_location"] = " | ".join(orf_trans_start_locations)
        
        #ORF_genome_start_location
        genome_location = result["ORF_genome_location"]
        strand = result["Strand"]
        orf_genome_start_locations = []
        for orf in genome_location.split(" | "):
            genome_matches = re.findall(r"\((\d+)-(\d+)\)", orf)
            if genome_matches:
                if strand == "+":
                    start_positions = [f"({m[0]})" for m in genome_matches]
                else:
                    #  When dealing with the negative strand, the interval is in the (end-start) format, so the start position should be taken from m[1]
                    start_positions = [f"({m[1]})" for m in genome_matches]
                orf_genome_start_locations.append(f"{orf.split('_')[0]}_[{'; '.join(start_positions)}]")
            else:
                orf_genome_start_locations.append(f"{orf.split('_')[0]}_[(NA)]")
        result["ORF_genome_start_location"] = " | ".join(orf_genome_start_locations)
    
    return results

def merge_original_with_results(original_data, results):
    """
    Merge the original data with the new results.
    """
    results_dict = {result["Trans_ID"]: result for result in results}
    
    original_df = pd.DataFrame(original_data)
    results_df = pd.DataFrame(results)
    
    merged_df = pd.merge(original_df, results_df, on="Trans_ID", how="left")
    
    merged_df.fillna({
        "Chromosome": "NA",
        "Strand": "NA",
        "ORF_genome_location": "NA",
        "ORF_trans_location": "NA",
        "ORF_genome_start_location": "NA",
        "ORF_trans_start_location": "NA"
    }, inplace=True)
    
    return merged_df

def write_output(output_file, merged_df):
    """
    Write the merged data into a CSV file, avoiding duplicate columns.
    """
    new_columns = [
        "Chromosome",
        "Strand",
        "ORF_genome_location",
        "ORF_trans_location",
        "ORF_genome_start_location",
        "ORF_trans_start_location"
    ]
    
    merged_df.to_csv(output_file, index=False)
    print(f"Intermediate results saved to {output_file}")

# -----------------------------
# Relative start position difference and mod3
# -----------------------------

def calculate_orf_difference(df):
    """
    Calculate the difference in relative start positions of ORFs and their modulo 3 remainder.
    Generate new columns: ORF_start_length_difference, ORF_trans_start_mod3.
    """
    def calculate_relative_start_position(orf_trans_start_location):
        if pd.isna(orf_trans_start_location) or orf_trans_start_location.strip() == "NA":
            return "", ""
        
        # Split each ORF entry by '|'
        orf_entries = [entry.strip() for entry in orf_trans_start_location.split('|')]
        formatted_entries = []
        mod3_entries = []
        
        for entry in orf_entries:
            # Extract ORF names (such as ORF5) and start positions
            orf_name = entry.split('_')[0]
            start_positions = re.findall(r'\((\d+)\)', entry)
            start_positions = list(map(int, start_positions))
            
            if start_positions:
                min_position = min(start_positions)
                # Calculate relative position differences (relative to the smallest start)
                relative_positions = [(pos - min_position + 1) for pos in start_positions]
                mod3_positions = [pos % 3 for pos in relative_positions]
                
                # Format the output
                formatted_entry = f"{orf_name}_[{'; '.join(f'({pos})' for pos in relative_positions)}]"
                mod3_entry = f"{orf_name}_[{'; '.join(f'({pos})' for pos in mod3_positions)}]"
                
                formatted_entries.append(formatted_entry)
                mod3_entries.append(mod3_entry)
        
        return '|'.join(formatted_entries), '|'.join(mod3_entries)
    
    results = df['ORF_trans_start_location'].apply(calculate_relative_start_position)
    df['ORF_start_length_difference'] = results.apply(lambda x: x[0])
    df['ORF_trans_start_mod3'] = results.apply(lambda x: x[1])
    
    return df

# -----------------------------
# Adjust genomic start positions + Calculate frame
# -----------------------------

def adjust_genome_start_and_frame(df):
    """
    Adjust genomic start positions and calculate frame information.
    Generate new columns: adj_ORF_genome_start_location, frame_info.
    """
    def process_orf_information(row):
        genome_start_str = row['ORF_genome_start_location']
        mod3_str = row['ORF_trans_start_mod3']
        chrom = row['Chromosome'] 

        if pd.isna(genome_start_str) or pd.isna(mod3_str) or genome_start_str.strip() == "NA" or mod3_str.strip() == "NA":
            return '', ''

        is_positive_strand = row['Strand'] == '+'

        orf_parts = [part.strip() for part in genome_start_str.split('|')]
        mod_parts = [part.strip() for part in mod3_str.split('|')]

        adjusted_orfs = []
        frame_infos = []

        for orf_part, mod_part in zip(orf_parts, mod_parts):
            orf_match = re.match(r'(ORF\d+)_\[(.+)\]', orf_part)
            if not orf_match:
                continue

            orf_name, positions = orf_match.groups()
            position_list = [int(pos) for pos in re.findall(r'\((\d+)\)', positions)]
            mod_list = [int(mod.strip("() ")) for mod in re.findall(r'\((\d+)\)', mod_part)]

            if len(position_list) != len(mod_list):
                continue

            adjusted_positions = []
            frame_values = []

            for start, mod in zip(position_list, mod_list):
                if is_positive_strand:
                    if mod == 1:
                        adjusted_start = start
                    elif mod == 2:
                        adjusted_start = start - 1
                    elif mod == 0:
                        adjusted_start = start - 2
                else:
                    if mod == 1:
                        adjusted_start = start
                    elif mod == 2:
                        adjusted_start = start + 1
                    elif mod == 0:
                        adjusted_start = start + 2

                adjusted_positions.append(f"({adjusted_start})")

                remainder = adjusted_start % 3
                if is_positive_strand:
                    if remainder == 1:
                        frame_values.append("(frame+1)")
                    elif remainder == 2:
                        frame_values.append("(frame+2)")
                    elif remainder == 0:
                        frame_values.append("(frame+3)")
                else:
                    if chrom in ['chr1', 'chr5', 'chr11', 'chr12', 'chr15', 'chr18', 'chr20', 'chrY']:
                        if remainder == 0:
                            frame_values.append("(frame-1)")
                        elif remainder == 2:
                            frame_values.append("(frame-2)")
                        elif remainder == 1:
                            frame_values.append("(frame-3)")
                    elif chrom in ['chr2', 'chr6', 'chr7','chr8', 'chr9', 'chr14' 'chr17', 'chr19', 'chr21', 'chr22']:
                        if remainder == 1:
                            frame_values.append("(frame-1)")
                        elif remainder == 0:
                            frame_values.append("(frame-2)")
                        elif remainder == 2:
                            frame_values.append("(frame-3)")
                    elif chrom in ['chr3', 'chr4', 'chr10', 'chr13', 'chr16', 'chrX']:
                        if remainder == 2:
                            frame_values.append("(frame-1)")
                        elif remainder == 1:
                            frame_values.append("(frame-2)")
                        elif remainder == 0:
                            frame_values.append("(frame-3)")
                    else:
                        if remainder == 0:
                            frame_values.append("(frame-1)")
                        elif remainder == 2:
                            frame_values.append("(frame-2)")
                        elif remainder == 1:
                            frame_values.append("(frame-3)")

            adjusted_orfs.append(f"{orf_name}_[{'; '.join(adjusted_positions)}]")
            frame_infos.append(f"{orf_name}_[{'; '.join(frame_values)}]")

        return '| '.join(adjusted_orfs), '| '.join(frame_infos)

    df[['adj_ORF_genome_start_location', 'frame_info']] = df.apply(
        lambda row: pd.Series(process_orf_information(row)),
        axis=1
    )

    return df


# -----------------------------
# Extract ORF genome end location, calculate differences and mod3, adjust end positions
# -----------------------------

def extract_orf_end_location(df):
    """
    Extract the ORF_genome_end_location column.
    Generate a new column: ORF_genome_end_location.
    """
    def extract_orf_end_location(row):
        if pd.isna(row['ORF_genome_location']) or row['ORF_genome_location'].strip() == "NA":
            # print(f"Skip rows with missing ORF_genome_location：{row.name}")
            return None
        
        orfs = [orf.strip() for orf in row['ORF_genome_location'].split('|')]
        strand = row['Strand']
        orf_end_locations = {}
    
        for orf in orfs:
            orf_name = orf.split('_')[0]
            intervals = orf.split(';')
            end_values = []
            
            for interval in intervals:
                match = re.search(r'\((\d+)-(\d+)\)', interval)
                if match:
                    start, end = int(match.group(1)), int(match.group(2))
                    end_value = max(start, end) if strand == '+' else min(start, end)
                    end_values.append(str(end_value))
            
            if orf_name in orf_end_locations:
                orf_end_locations[orf_name].extend(end_values)
            else:
                orf_end_locations[orf_name] = end_values
    
        formatted_end_locations = [f"{orf}_[{'; '.join([f'({end})' for end in ends])}]" for orf, ends in orf_end_locations.items()]
        return '|'.join(formatted_end_locations).replace('| ', '|')
    
    df['ORF_genome_end_location'] = df.apply(extract_orf_end_location, axis=1)
    
    return df

def calculate_orf_end_length_difference(df):
    """
    Calculate the ORF_end_length_difference column.
    """
    def calculate_orf_end_length_difference(row):
        if pd.isna(row['ORF_genome_end_location']) or pd.isna(row['adj_ORF_genome_start_location']) \
           or row['ORF_genome_end_location'].strip() == "NA" or row['adj_ORF_genome_start_location'].strip() == "NA":
            return None
        
        orf_end_groups = [group.strip() for group in row['ORF_genome_end_location'].split('|')]
        start_groups = [group.strip() for group in row['adj_ORF_genome_start_location'].split('|')]
        strand = row['Strand']
        orf_differences = {}
    
        for orf_end_group, start_group in zip(orf_end_groups, start_groups):
            orf_name = orf_end_group.split('_')[0]
            orf_ends = re.findall(r'\((\d+)\)', orf_end_group)
            starts = re.findall(r'\((\d+)\)', start_group)
            
            if len(orf_ends) != len(starts):
                print(f"Row {row.name} ORF {orf_name} position counts do not match")
                return None
    
            differences = []
            for orf_end, adj_start in zip(orf_ends, starts):
                orf_end, adj_start = int(orf_end), int(adj_start)
                if strand == '+':
                    difference = orf_end - adj_start + 1
                else:
                    difference = adj_start - orf_end + 1
                differences.append(str(difference))
    
            if orf_name in orf_differences:
                orf_differences[orf_name].extend(differences)
            else:
                orf_differences[orf_name] = differences
    
        formatted_differences = [f"{orf}_[{'; '.join([f'({diff})' for diff in diffs])}]" for orf, diffs in orf_differences.items()]
        return '|'.join(formatted_differences).replace('| ', '|')
    
    df['ORF_end_length_difference'] = df.apply(calculate_orf_end_length_difference, axis=1)
    
    return df

def calculate_orf_trans_end_mod3(df):
    """
    Calculate ORF_trans_end_mod3 column.
    """
    def calculate_orf_trans_end_mod3(row):
        if pd.isna(row['ORF_end_length_difference']) or row['ORF_end_length_difference'].strip() == "NA":
            return None
        
        orf_diff_groups = [group.strip() for group in row['ORF_end_length_difference'].split('|')]
        mod3_results = {}
    
        for orf_diff_group in orf_diff_groups:
            orf_name = orf_diff_group.split('_')[0]
            lengths = re.findall(r'\((\d+)\)', orf_diff_group)
            mods = [str(int(length) % 3) for length in lengths]
    
            if orf_name in mod3_results:
                mod3_results[orf_name].extend(mods)
            else:
                mod3_results[orf_name] = mods
    
        formatted_mod3 = [f"{orf}_[{'; '.join([f'({mod})' for mod in mods])}]" for orf, mods in mod3_results.items()]
        return '|'.join(formatted_mod3).replace('| ', '|')
    
    df['ORF_trans_end_mod3'] = df.apply(calculate_orf_trans_end_mod3, axis=1)
    
    return df

def calculate_adj_orf_genome_end_location(df):
    """
    Calculate adj_ORF_genome_end_location column.
    """
    def calculate_adj_orf_genome_end_location(row):
        if pd.isna(row['ORF_genome_end_location']) or pd.isna(row['ORF_trans_end_mod3']) \
           or row['ORF_genome_end_location'].strip() == "NA" or row['ORF_trans_end_mod3'].strip() == "NA":
            return None
        
        orf_end_groups = [group.strip() for group in row['ORF_genome_end_location'].split('|')]
        mod3_groups = [group.strip() for group in row['ORF_trans_end_mod3'].split('|')]
        strand = row['Strand']
        adj_end_locations = {}
    
        for orf_end_group, mod3_group in zip(orf_end_groups, mod3_groups):
            orf_name = orf_end_group.split('_')[0]
            orf_ends = re.findall(r'\((\d+)\)', orf_end_group)
            mods = re.findall(r'\((\d+)\)', mod3_group)
            
            if len(orf_ends) != len(mods):
                print(f"Row {row.name} ORF {orf_name} position counts do not match")
                return None
    
            adjusted_ends = []
            for orf_end, mod in zip(orf_ends, mods):
                orf_end, mod = int(orf_end), int(mod)
                if strand == '+':
                    adj_end = orf_end - mod
                else:
                    adj_end = orf_end + mod
                adjusted_ends.append(str(adj_end))
    
            if orf_name in adj_end_locations:
                adj_end_locations[orf_name].extend(adjusted_ends)
            else:
                adj_end_locations[orf_name] = adjusted_ends
    
        formatted_adj_end = [f"{orf}_[{'; '.join([f'({adj_end})' for adj_end in ends])}]" for orf, ends in adj_end_locations.items()]
        return '|'.join(formatted_adj_end).replace('| ', '|')
    
    df['adj_ORF_genome_end_location'] = df.apply(calculate_adj_orf_genome_end_location, axis=1)
    
    return df

def calculate_adj_orf_genome_location(df):
    """
    Calculate adj_ORF_genome_location column.
    """
    def calculate_adj_orf_genome_location(row):
        if pd.isna(row['adj_ORF_genome_start_location']) or pd.isna(row['adj_ORF_genome_end_location']) \
           or row['adj_ORF_genome_start_location'].strip() == "NA" or row['adj_ORF_genome_end_location'].strip() == "NA":
            return None
        
        orf_start_groups = [group.strip() for group in row['adj_ORF_genome_start_location'].split('|')]
        orf_end_groups = [group.strip() for group in row['adj_ORF_genome_end_location'].split('|')]
        strand = row['Strand']
        adj_locations = {}
    
        for orf_start_group, orf_end_group in zip(orf_start_groups, orf_end_groups):
            orf_name = orf_start_group.split('_')[0]
            starts = re.findall(r'\((\d+)\)', orf_start_group)
            ends = re.findall(r'\((\d+)\)', orf_end_group)
            
            if len(starts) != len(ends):
                print(f"Row {row.name} ORF {orf_name} position counts do not match")
                return None
    
            location_pairs = [f"({start}-{end})" for start, end in zip(starts, ends)]
            if orf_name in adj_locations:
                adj_locations[orf_name].extend(location_pairs)
            else:
                adj_locations[orf_name] = location_pairs
    
        formatted_adj_location = [f"{orf}_[{'; '.join(locations)}]" for orf, locations in adj_locations.items()]
        return '|'.join(formatted_adj_location).replace('| ', '|')
    
    df['adj_ORF_genome_location'] = df.apply(calculate_adj_orf_genome_location, axis=1)
    
    return df

# -----------------------------
# Determine WIG file paths and convert to bedGraph format
# -----------------------------

def determine_wig_files_with_format(row, base_directory):
    """
    Based on frame_info, determine the WIG file paths.
    Generate a new column: phylocsf_file.
    """
    if pd.isna(row['Chromosome']) or pd.isna(row['frame_info']) or row['frame_info'].strip() == "NA":
        return ""  # Return an empty string if data is missing
    
    chromosome_folder = f"output_{row['Chromosome'].strip()}/"  # Folder based on chromosome
    frame_info = str(row['frame_info'])  # Convert to string to handle NaN values
    orfs = frame_info.split('|')
    orf_files = []

    # Process each ORF
    for orf_entry in orfs:
        orf_entry = orf_entry.strip()
        orf_match = re.match(r"(ORF\d+)_\[(.*)\]", orf_entry)
        if orf_match:
            orf_id = orf_match.group(1)
            frames = orf_match.group(2).split(';')
            frame_paths = []
            for frame in frames:
                frame = frame.strip("() ")
                frame_match = re.match(r"frame([+-]\d)", frame)
                if frame_match:
                    frame_suffix = frame_match.group(1)
                    wig_file_path = os.path.join(base_directory, chromosome_folder, f"PhyloCSFRaw{frame_suffix}.bw")
                    frame_paths.append(f"({wig_file_path})")
                else:
                    print(f"  Frame pattern not matched: {frame}")
            orf_files.append(f"{orf_id}_[{'; '.join(frame_paths)}]")
        else:
            print(f"ORF pattern not matched: {orf_entry}")
    return '|'.join(orf_files)

def convert_to_bedgraph_format(row):
    """
    Convert adj_ORF_genome_location from 1-based to 0-based bedGraph format. 
    Generate a new column: adj_ORF_genome_location_bedgraph.
    """
    if pd.isna(row['adj_ORF_genome_location']) or row['adj_ORF_genome_location'].strip() == "NA":
        return ""  # Return an empty string if data is missing

    adj_location = str(row['adj_ORF_genome_location'])
    orfs = adj_location.split('|')
    bedgraph_entries = []

    # Process each ORF entry
    for orf_entry in orfs:
        orf_entry = orf_entry.strip()
        orf_match = re.match(r"(ORF\d+)_\[(.*)\]", orf_entry)
        if orf_match:
            orf_id = orf_match.group(1)
            locations = orf_match.group(2).split(';')
            location_entries = []
            for loc in locations:
                loc = loc.strip("() ")
                start_end_match = re.match(r"(\d+)-(\d+)", loc)
                if start_end_match:
                    start = int(start_end_match.group(1)) - 1  # Convert to 0-based by subtracting 1
                    end = int(start_end_match.group(2))
                    location_entries.append(f"({start}-{end})")
                else:
                    print(f"  Location pattern not matched: {loc}")
            bedgraph_entries.append(f"{orf_id}_[{'; '.join(location_entries)}]")
        else:
            print(f"ORF pattern not matched: {orf_entry}")
    return '|'.join(bedgraph_entries).replace('| ', '|')

def determine_wig_files(df, base_directory):
    """
    Determine WIG file paths and generate the phylocsf_file column.
    """
    df['phylocsf_file'] = df.apply(determine_wig_files_with_format, axis=1, base_directory=base_directory)
    return df

def convert_adj_genome_location_to_bedgraph(df):
    """
    Convert adj_ORF_genome_location to bedGraph format. 
    Generate a new column: adj_ORF_genome_location_bedgraph
    """
    df['adj_ORF_genome_location_bedgraph'] = df.apply(convert_to_bedgraph_format, axis=1)
    return df

# -----------------------------
# Query BigWig files to obtain PhyloCSF scores
# -----------------------------

def query_bigwig(bigwig_file, chrom, start, end):
    """
    Use the bigWigSummary tool to query BigWig files.
    
    Parameters:
    - bigwig_file (str): The path to the BigWig file.
    - chrom (str): The chromosome name.
    - start (int): The start position.
    - end (int): The end position.
    
    Returns:
    - float or None: The PhyloCSF score retrieved by the query, returns None if the query fails.
    """
    if start > end:
        start, end = end, start

    print(f"Querying: {bigwig_file}, {chrom}:{start}-{end}")  # Debug: Print each query
    try:
        result = subprocess.check_output([
            "bigWigSummary", bigwig_file, chrom, str(start), str(end), "1"
        ])
        return float(result.strip())  # Convert the result to float and strip newline
    except subprocess.CalledProcessError as e:
        print(f"Error querying {bigwig_file} for {chrom}:{start}-{end}: {e}")
        return None
    except ValueError:
        print(f"No data for {chrom}:{start}-{end} in {bigwig_file}")
        return None

def get_phylocsf_scores(row):
    """
    Process each ORF and retrieve PhyloCSF scores.
    Generate new columns: ORF_phylocsf_score, ORF_phylocsf_merge_score
    """
    if pd.isna(row['phylocsf_file']) or pd.isna(row['adj_ORF_genome_location_bedgraph']) or pd.isna(row['Chromosome']):
        print("Skipping row due to missing values")  # Debug: Missing data skip
        return "", []
    
    phylocsf_entries = row['phylocsf_file'].split('|')
    bedgraph_entries = row['adj_ORF_genome_location_bedgraph'].split('|')
    orf_scores = []

    for phylocsf_entry, bedgraph_entry in zip(phylocsf_entries, bedgraph_entries):
        orf_match = re.match(r"(ORF\d+)_\[(.*)\]", bedgraph_entry.strip())
        file_match = re.match(r"(ORF\d+)_\[(.*)\]", phylocsf_entry.strip())

        if orf_match and file_match:
            orf_id = orf_match.group(1)
            locations = orf_match.group(2).split(';')
            files = file_match.group(2).split(';')
            location_scores = []

            for loc, bw_file in zip(locations, files):
                loc = loc.strip("() ").strip()
                bw_file = bw_file.strip("() ").strip()

                chrom = row['Chromosome']
                start_end_match = re.match(r"(\d+)-(\d+)", loc)
                if start_end_match:
                    start = int(start_end_match.group(1))
                    end = int(start_end_match.group(2))
                    print(f"Processing ORF {orf_id} on {chrom} from {start} to {end} in {bw_file}")  # Debug

                    score = query_bigwig(bw_file, chrom, start, end)
                    if score is not None:
                        location_scores.append(score)
                    else:
                        location_scores.append(None)

            # Format the scores as required in ORF_phylocsf_score format
            formatted_scores = "; ".join(f"({s})" if s is not None else "(NA)" for s in location_scores)
            orf_scores.append(f"{orf_id}_[{formatted_scores}]")
        else:
            print(f"ORF or file pattern not matched: {phylocsf_entry}, {bedgraph_entry}")

    return '|'.join(orf_scores), orf_scores

def calculate_merged_scores(orf_scores):
    """
    Calculate and format the merged score, ignoring NA values.
    
    Parameters:
    - orf_scores (list): A list of ORF scores.
    
    Returns:
    - str: The formatted merged score string.
    """
    merged_scores = []

    for orf_score in orf_scores:
        orf_match = re.match(r"(ORF\d+)_\[(.*)\]", orf_score.strip())
        if orf_match:
            orf_id = orf_match.group(1)
            # Use regex to extract numbers inside parentheses
            values = re.findall(r"\((-?\d+\.?\d*)\)", orf_match.group(2))
            try:
                # Convert extracted values to floats and filter out NA values
                scores = [float(value) for value in values if value != "NA"]
                avg_score = mean(scores) if scores else None
                formatted_avg = f"{orf_id}_[{avg_score if avg_score is not None else 'NA'}]"
                merged_scores.append(formatted_avg)
            except ValueError as e:
                print(f"Error converting score to float in {orf_score}: {e}")  # Debug output
                merged_scores.append(f"{orf_id}_[NA]")  # Add NA if conversion fails
        else:
            print(f"Invalid ORF score format: {orf_score}")

    return '|'.join(merged_scores)

def check_positive_orf_scores(orf_merged_scores):
    """
    Check if all ORF scores are positive, mark as UNKNOWN if there is any NA.
    
    Parameters:
    - orf_merged_scores (str): The merged ORF score string.
    
    Returns:
    - str: "TRUE", "FALSE" or "UNKNOWN".
    """
    if pd.isna(orf_merged_scores):
        return "UNKNOWN"
    orf_entries = orf_merged_scores.split('|')
    for orf_entry in orf_entries:
        scores = re.findall(r"\[(-?\d+\.?\d*)\]", orf_entry)
        if "NA" in orf_entry:
            return "UNKNOWN"
        if any(float(score) <= 0 for score in scores):
            return "FALSE"
    return "TRUE"

def add_phylocsf_scores(df):
    """
    Add new columns to your DataFrame for storing PhyloCSF scores.
    Populate these columns with the relevant PhyloCSF scores for each ORF.
    """
    print("Processing PhyloCSF scores...")
    results = df.apply(get_phylocsf_scores, axis=1, result_type="expand")
    df['ORF_phylocsf_score'] = results[0]
    df['ORF_phylocsf_merge_score'] = results[1].apply(calculate_merged_scores)

    print("Checking if all ORF scores are positive...")
    df['ORF_all_positive'] = df['ORF_phylocsf_merge_score'].apply(check_positive_orf_scores)

    return df

# -----------------------------
# Main Function
# -----------------------------

def main(gff_file, orf_csv_file, output_file, base_directory):
    """
    Main Function：
    1. Parse the GFF3 file.
    2. Parse the ORF CSV file.
    3. Map ORFs to genomic coordinates.
    4. Extract ORF start position information.
    5. Merge data.
    6. Write intermediate results.
    7. Use Pandas for further calculations.
    8. Determine WIG file paths and convert to bedGraph format.
    9. Query BigWig files to obtain PhyloCSF scores.
    10. Calculate and format the merged PhyloCSF scores.
    11. Check if all ORF scores are positive.
    12. Filter and save the final results.
    """
    print("Step 1/12: Parsing GFF3 file...")
    transcript_exons = parse_gff(gff_file)
    print(f"Parsed {len(transcript_exons)} transcripts with exon information.")
    
    print("Step 2/12: Parsing ORF CSV file...")
    orf_data, original_data = parse_orf_csv(orf_csv_file)
    print(f"Parsed {len(orf_data)} transcripts with ORF data.")
    
    print("Step 3/12: Mapping ORFs to genome coordinates...")
    results = map_orf_to_genome(transcript_exons, orf_data)
    print(f"Mapped {len(results)} transcripts to genome coordinates.")
    
    print("Step 4/12: Extracting ORF start locations...")
    results = extract_orf_start_locations(results)
    
    print("Step 5/12: Merging original data with new results...")
    merged_df = merge_original_with_results(original_data, results)
    
    print("Step 6/12: Writing intermediate CSV...")
    write_output(output_file, merged_df)
    
    # Read intermediate results into a DataFrame for further calculations.
    print("Step 7/12: Loading intermediate CSV into DataFrame for further calculations...")
    try:
        df = pd.read_csv(output_file)
        print(f"Successfully loaded intermediate file: {output_file}")
    except Exception as e:
        print(f"Failed to load intermediate file: {e}")
        sys.exit(1)
    
    print("Step 8/12: Calculating ORF start length differences and mod3 info...")
    df = calculate_orf_difference(df)
    
    print("Step 9/12: Adjusting genome start locations and calculating frame info...")
    df = adjust_genome_start_and_frame(df)
    
    print("Step 10/12: Extracting ORF genome end locations...")
    df = extract_orf_end_location(df)
    
    print("Step 11/12: Calculating ORF end length differences and mod3 info...")
    df = calculate_orf_end_length_difference(df)
    df = calculate_orf_trans_end_mod3(df)
    
    print("Step 12/12: Calculating adjusted genome end locations...")
    df = calculate_adj_orf_genome_end_location(df)
    
    print("Step 13/12: Calculating adjusted genome location...")
    df = calculate_adj_orf_genome_location(df)
    
    print("Step 14/12: Determining wig file paths and converting to bedgraph format...")
    df = determine_wig_files(df, base_directory)
    df = convert_adj_genome_location_to_bedgraph(df)
    
    print("Step 15/12: Adding PhyloCSF scores...")
    df = add_phylocsf_scores(df)
    
    # Filter out rows where phylocsf_file is not null.
    print("Step 16/12: Filtering rows with non-empty phylocsf_file...")
    df_filtered = df[df['phylocsf_file'].str.len() > 0]
    
    # Output the first few rows of the new columns.
    print("Generated new columns (first few rows):")
    print(df_filtered[['ORF_genome_end_location', 'ORF_end_length_difference', 'ORF_trans_end_mod3', 
                      'adj_ORF_genome_end_location', 'adj_ORF_genome_location', 'frame_info', 
                      'phylocsf_file', 'adj_ORF_genome_location_bedgraph',
                      'ORF_phylocsf_score', 'ORF_phylocsf_merge_score', 'ORF_all_positive']].head())
    
    # Save the final results to a specified output file.
    print("Saving final results to output file...")
    try:
        df_filtered.to_csv(output_file, index=False)
        print(f"Processing complete. Final results saved to '{output_file}'.")
    except Exception as e:
        print(f"Failed to save final output: {e}")

# -----------------------------
# Command-line Interface
# -----------------------------
if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(
        description=(
            "1) Map ORF transcript positions to genome coordinates using GFF3 and ORF CSV files, "
            "2) Calculate ORF start-length differences and mod3 remainder, "
            "3) Adjust genome start positions and compute frame info, "
            "4) Extract ORF genome end locations, "
            "5) Calculate ORF end-length differences and mod3 remainder, "
            "6) Adjust genome end positions and compute final genome coordinates, "
            "7) Determine wig file paths and convert to bedgraph format, "
            "8) Query BigWig files to retrieve PhyloCSF scores, "
            "9) Calculate and merge PhyloCSF scores, "
            "10) Check if all ORF scores are positive."
        )
    )
    parser.add_argument('--gff3_file', required=True, help="Path to GFF3 file")
    parser.add_argument('--orf_csv_file', required=True, help="Path to ORF CSV file")
    parser.add_argument('--output_file', required=True, help="Path to final output CSV file")
    parser.add_argument('--base_directory', required=True, help="Base directory for chromosome-specific output folders")
    args = parser.parse_args()
    
    main(args.gff3_file, args.orf_csv_file, args.output_file, args.base_directory)
