# -*- coding: utf-8 -*-
__author__ = "legendzdy@dingtalk.com"
"""
Author: legendzdy@dingtalk.com
Date: 20250316
Description: This file contains utilities for gene annotation from GTF files and isoform ID parsing.
"""
from typing import Dict, List, Tuple
import pandas as pd

def parse_gtf(gtf_file: str) -> tuple[Dict[str, str], Dict[str, str], Dict[str, tuple]]:
    """
    Parse GTF file to create mapping from gene_id to gene_name and transcript_id to gene_id.
    Returns tuple (gene_map, transcript_map, gene_coords) where:
    - gene_map: {gene_id: gene_name}
    - transcript_map: {transcript_id: gene_id}
    - gene_coords: {gene_id: (chrom, start, end, strand)} with chrom without 'chr' prefix
    """
    gene_map = {}
    transcript_map = {}
    gene_coords = {}
    with open(gtf_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            if len(fields) < 9:
                continue
            chrom = fields[0]
            # Remove 'chr' prefix if present for consistent matching
            if chrom.lower().startswith('chr'):
                chrom = chrom[3:]
            feature_type = fields[2]
            start = int(fields[3])
            end = int(fields[4])
            strand = fields[6]
            attributes = fields[8]
            # Parse key-value pairs
            gene_id = None
            gene_name = None
            transcript_id = None
            parts = attributes.split(';')
            for part in parts:
                part = part.strip()
                if not part:
                    continue
                if ' ' in part:
                    key, value = part.split(' ', 1)
                    key = key.strip()
                    value = value.strip().strip('"')
                    if key == 'gene_id':
                        gene_id = value
                    elif key == 'gene_name':
                        gene_name = value
                    elif key == 'transcript_id':
                        transcript_id = value
            if gene_id:
                gene_map[gene_id] = gene_name if gene_name else gene_id
                if transcript_id:
                    transcript_map[transcript_id] = gene_id
                # Store coordinates for gene features only
                if feature_type == 'gene':
                    gene_coords[gene_id] = (chrom, start, end, strand)
    return gene_map, transcript_map, gene_coords


def extract_gene_ids_from_isoform_str(isoform_str: str, gene_map: Dict[str, str], transcript_map: Dict[str, str]) -> Tuple[str, str]:
    """
    Extract gene IDs and names from isoform_ids string.
    
    Args:
        isoform_str: comma-separated string of isoform_id_gene_id pairs
        gene_map: dict mapping gene_id -> gene_name
        transcript_map: dict mapping transcript_id -> gene_id
    
    Returns:
        tuple: (gene_id_str, gene_name_str)
        - gene_id_str: ENS gene IDs joined by '|' (empty if none found)
        - gene_name_str: corresponding gene names joined by '|'
    """
    # Split by comma, then each element split by underscore
    isoforms = [item.strip() for item in isoform_str.split(',')]
    # Collect ENS gene IDs from all isoforms in this row
    ens_gene_ids = []
    for item in isoforms:
        parts = item.split('_')
        if len(parts) >= 2:
            # Format: uuid_geneId or transcriptId_geneId
            gene_id = parts[-1]  # gene_id is after last underscore
            # Accept gene ID if it's in gene_map or looks like an ENSEMBL ID
            if gene_id in gene_map or gene_id.startswith('ENS'):
                ens_gene_ids.append(gene_id)
        else:
                # No underscore: could be transcript ID, gene ID, or UUID
                if item.startswith('ENS'):
                    # Generic handling for any species ENSEMBL IDs
                    # First check if it's a transcript ID (mapped to gene)
                    gene_id = transcript_map.get(item)
                    if gene_id:
                        # This is a transcript ID that maps to a gene
                        ens_gene_ids.append(gene_id)
                    elif item in gene_map:
                        # Direct gene ID
                        ens_gene_ids.append(item)
                    else:
                        # ID not found in maps, try to infer type from pattern
                        # ENSEMBL IDs format: ENS[spp][type][number] where type is G (gene), T (transcript), etc.
                        # Find the first occurrence of G or T after ENS prefix
                        # This is a heuristic and may not be perfect for all cases
                        if 'G' in item and ('T' not in item or item.find('G') < item.find('T')):
                            # Likely a gene ID (contains G and either no T or G comes before T)
                            # But we can't map it since it's not in gene_map, so ignore
                            pass
                        # Other ENS prefixes (e.g., ENSMUSP, ENSE, etc.) ignore
    # Remove duplicates while preserving order
    unique_ens_gene_ids = []
    seen = set()
    for g in ens_gene_ids:
        if g not in seen:
            seen.add(g)
            unique_ens_gene_ids.append(g)
    
    # Join gene IDs with | if multiple, else single or empty
    if unique_ens_gene_ids:
        gene_id_str = '|'.join(unique_ens_gene_ids)
        # Get corresponding gene names in same order
        name_parts = []
        for g in unique_ens_gene_ids:
            name = gene_map.get(g, g)  # Use gene_id if name not found
            name_parts.append(name)
        gene_name_str = '|'.join(name_parts)
    else:
        gene_id_str = ''
        gene_name_str = ''
    
    return gene_id_str, gene_name_str


def parse_feature_id(feature_id: str) -> tuple[str, int, int]:
    """
    Parse feature_id string to extract chromosome, start, end coordinates.
    Handles formats:
      inclusion_CHR:start-end
      inclusion_CHR:start
      inclusion_CHR:start-end-suffix (suffix ignored)
    Returns (chrom, start, end) where chrom does not have 'chr' prefix.
    """
    # Remove inclusion_ or exclusion_ prefix
    if feature_id.startswith('inclusion_'):
        coord_str = feature_id[len('inclusion_'):]
    elif feature_id.startswith('exclusion_'):
        coord_str = feature_id[len('exclusion_'):]
    else:
        # If no prefix, assume entire string is coordinate
        coord_str = feature_id
    
    # Split chromosome and position part
    if ':' not in coord_str:
        raise ValueError(f"Invalid feature_id format: {feature_id}")
    chrom, pos_part = coord_str.split(':', 1)
    # Remove 'chr' prefix if present
    if chrom.lower().startswith('chr'):
        chrom = chrom[3:]
    
    # Parse position part which may contain range and/or suffix
    # Split by '-'
    parts = pos_part.split('-')
    # First part is start position
    try:
        start = int(parts[0])
    except ValueError:
        raise ValueError(f"Invalid start position in feature_id: {feature_id}")
    
    # Determine end position
    if len(parts) >= 2:
        # Check if second part is numeric (could be end position or suffix)
        if parts[1].isdigit():
            candidate = int(parts[1])
            # If candidate is greater than start, assume it's end position; otherwise treat as suffix
            if candidate > start:
                end = candidate
            else:
                # Suffix (e.g., -1), treat as single position event
                end = start
        else:
            # Suffix, not a position, treat as single position event
            end = start
    else:
        end = start
    
    return chrom, start, end


def annotate_by_coordinates(input_df: pd.DataFrame, gene_coords: Dict[str, tuple], gene_map: Dict[str, str]) -> pd.DataFrame:
    """
    Annotate alternative splicing events using genomic coordinates from feature_id.
    
    Args:
        input_df: pandas DataFrame with 'feature_id' column
        gene_coords: dict mapping gene_id -> (chrom, start, end, strand)
        gene_map: dict mapping gene_id -> gene_name
    
    Returns:
        pandas DataFrame with added 'gene_id' and 'gene_name' columns
    """
    if 'feature_id' not in input_df.columns:
        raise ValueError("Input DataFrame must contain 'feature_id' column")
    
    # Build chromosome-indexed gene intervals
    chrom_genes = {}
    for gene_id, (chrom, start, end, strand) in gene_coords.items():
        chrom_genes.setdefault(chrom, []).append((start, end, gene_id))
    
    gene_ids = []
    gene_names = []
    
    for idx, feature_id in enumerate(input_df['feature_id']):
        try:
            chrom, start, end = parse_feature_id(feature_id)
        except ValueError as e:
            # If parsing fails, leave empty
            gene_ids.append('')
            gene_names.append('')
            continue
        
        # Find overlapping genes on this chromosome
        overlapping_gene_ids = []
        if chrom in chrom_genes:
            for gene_start, gene_end, gene_id in chrom_genes[chrom]:
                if gene_start <= end and gene_end >= start:
                    overlapping_gene_ids.append(gene_id)
        
        # Remove duplicates while preserving order
        unique_gene_ids = []
        seen = set()
        for g in overlapping_gene_ids:
            if g not in seen:
                seen.add(g)
                unique_gene_ids.append(g)
        
        if unique_gene_ids:
            gene_id_str = '|'.join(unique_gene_ids)
            # Get gene names
            name_parts = []
            for g in unique_gene_ids:
                name = gene_map.get(g, g)
                name_parts.append(name)
            gene_name_str = '|'.join(name_parts)
        else:
            gene_id_str = ''
            gene_name_str = ''
        
        gene_ids.append(gene_id_str)
        gene_names.append(gene_name_str)
    
    # Insert new columns after feature_id column if exists, otherwise at end
    if 'feature_id' in input_df.columns:
        feature_idx = list(input_df.columns).index('feature_id')
        input_df.insert(feature_idx + 1, 'gene_id', gene_ids)  # type: ignore
        input_df.insert(feature_idx + 2, 'gene_name', gene_names)  # type: ignore
    else:
        input_df['gene_id'] = gene_ids
        input_df['gene_name'] = gene_names
    
    return input_df


def annotate_isoform_events(input_df: pd.DataFrame, gene_map: Dict[str, str], transcript_map: Dict[str, str]) -> pd.DataFrame:
    """
    Annotate isoform events dataframe with gene IDs and names.
    
    Args:
        input_df: pandas DataFrame with 'isoform_ids' column
        gene_map: dict mapping gene_id -> gene_name
        transcript_map: dict mapping transcript_id -> gene_id
    
    Returns:
        pandas DataFrame with added 'gene_id' and 'gene_name' columns
    """
    if 'isoform_ids' not in input_df.columns:
        raise ValueError("Input DataFrame must contain 'isoform_ids' column")
    
    gene_ids = []
    gene_names = []
    warning_count = 0
    
    for idx, isoform_str in enumerate(input_df['isoform_ids']):
        gene_id_str, gene_name_str = extract_gene_ids_from_isoform_str(isoform_str, gene_map, transcript_map)
        
        # Check for multiple distinct gene IDs
        if gene_id_str and '|' in gene_id_str:
            if warning_count < 5:  # Limit warning messages
                print(f"Warning: row {idx+1} has multiple distinct ENS gene IDs: {gene_id_str.split('|')}. All will be kept.")
            warning_count += 1
        
        gene_ids.append(gene_id_str)
        gene_names.append(gene_name_str)
    
    if warning_count > 0:
        print(f"Total rows with multiple ENS gene IDs: {warning_count}")
    
    # Insert new columns after isoform_ids column
    isoform_idx = list(input_df.columns).index('isoform_ids')
    input_df.insert(isoform_idx + 1, 'gene_id', gene_ids)  # type: ignore
    input_df.insert(isoform_idx + 2, 'gene_name', gene_names)  # type: ignore
    
    return input_df


def annotate_single_file(input_file: str, gene_map: dict, transcript_map: dict, gene_coords: dict, output_file: str = "", progress=None) -> str:
    """
    Annotate a single TSV file and return output path.
    Uses coordinate-based annotation if feature_id column exists, otherwise isoform-based annotation.
    """
    from pathlib import Path
    import pandas as pd
    import os

    if not output_file:
        input_path = Path(input_file)
        output_file = str(input_path.with_suffix('.annotated.tsv'))

    if progress:
        progress.add_task(description=f"Processing {os.path.basename(input_file)}...", total=None)

    # Ensure output directory exists
    output_path = Path(output_file)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(input_file, sep='\t')
    
    # Choose annotation method based on available columns
    if 'feature_id' in df.columns:
        df = annotate_by_coordinates(df, gene_coords, gene_map)
    elif 'isoform_ids' in df.columns:
        df = annotate_isoform_events(df, gene_map, transcript_map)
    else:
        raise ValueError("Input file must contain either 'feature_id' or 'isoform_ids' column")
    
    df.to_csv(output_file, sep='\t', index=False)

    return output_file


def count_inclusion_exclusion(df: pd.DataFrame) -> tuple[int, int]:
    """
    Count inclusion and exclusion events based on feature_id prefix.
    Returns (inclusion_count, exclusion_count)
    """
    if 'feature_id' not in df.columns:
        return 0, 0
    inclusion = df['feature_id'].str.startswith('inclusion_').sum()
    exclusion = df['feature_id'].str.startswith('exclusion_').sum()
    return inclusion, exclusion


def plot_splicing_counts(counts_dict: dict, output_prefix: str):
    """
    Generate stacked bar chart for inclusion/exclusion counts across splicing types.
    Horizontal stacked bar chart with Y-axis as splicing types and X-axis as counts.
    counts_dict: {splicing_type: (inclusion_count, exclusion_count)}
    """
    import sys
    try:
        import matplotlib.pyplot as plt
        import numpy as np
    except ImportError:
        print("Warning: matplotlib not installed. Skipping plot generation.", file=sys.stderr)
        return
    
    type_order = ['IR', 'Alt3', 'Alt5', 'ES']
    types = []
    for t in type_order:
        if t in counts_dict:
            types.append(t)
    # Add any remaining types not in predefined order
    for t in counts_dict.keys():
        if t not in types:
            types.append(t)
    
    inclusion_counts = [counts_dict[t][0] for t in types]
    exclusion_counts = [counts_dict[t][1] for t in types]
    
    y = np.arange(len(types))
    height = 0.6
    
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.barh(y, inclusion_counts, height, label='Inclusion', color='#DD95B9')
    ax.barh(y, exclusion_counts, height, left=inclusion_counts, label='Exclusion', color='#6A90CA')
    
    ax.set_xlabel('Count', fontsize=14)
    ax.set_ylabel('Splicing Type', fontsize=14)
    ax.set_title('Inclusion/Exclusion Events by Splicing Type', fontsize=16)
    ax.set_yticks(y)
    ax.set_yticklabels(types)
    ax.legend()
    
    # Add value labels on bars
    for i, (inc, exc) in enumerate(zip(inclusion_counts, exclusion_counts)):
        total = inc + exc
        if total > 0:
            # Inclusion label
            if inc > 0:
                ax.text(inc/2, i, str(inc), ha='center', va='center', color='white', fontweight='bold')
            # Exclusion label
            if exc > 0:
                ax.text(inc + exc/2, i, str(exc), ha='center', va='center', color='white', fontweight='bold')
    
    plt.tight_layout()
    
    # Ensure output directory exists
    from pathlib import Path
    plot_dir = Path(output_prefix).parent
    plot_dir.mkdir(parents=True, exist_ok=True)
    
    base_path = f"{output_prefix}_splicing_counts"
    saved_files = []
    
    # Save PNG (default format)
    try:
        png_path = f"{base_path}.png"
        plt.savefig(png_path, dpi=300)
        saved_files.append(png_path)
    except Exception as e:
        print(f"Warning: Could not save PNG format: {e}", file=sys.stderr)
    
    # Save PDF
    try:
        pdf_path = f"{base_path}.pdf"
        plt.savefig(pdf_path, dpi=300)
        saved_files.append(pdf_path)
    except Exception as e:
        print(f"Warning: Could not save PDF format: {e}", file=sys.stderr)
    
    # Try to save TIFF (requires pillow)
    try:
        tiff_path = f"{base_path}.tiff"
        plt.savefig(tiff_path, dpi=300)
        saved_files.append(tiff_path)
    except Exception as e:
        print(f"Warning: Could not save TIFF format: {e}", file=sys.stderr)
        print("Tip: Install pillow package for TIFF support: pip install pillow", file=sys.stderr)
    
    plt.close()
    
    if saved_files:
        print(f"Plots saved to:")
        for file_path in saved_files:
            print(f"  {file_path}")
    else:
        print("Warning: No plot formats were successfully saved.", file=sys.stderr)