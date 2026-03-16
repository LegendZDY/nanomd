# -*- coding: utf-8 -*-
__author__ = "legendzdy@dingtalk.com"
"""
Author: legendzdy@dingtalk.com
Date: 20250316
Description: This file contains utilities for gene annotation from GTF files and isoform ID parsing.
"""
from typing import Dict, List, Tuple
import pandas as pd

def parse_gtf(gtf_file: str) -> Dict[str, str]:
    """
    Parse GTF file to create mapping from gene_id to gene_name.
    Returns dict {gene_id: gene_name}
    """
    gene_map = {}
    with open(gtf_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            if len(fields) < 9:
                continue
            attributes = fields[8]
            # Parse key-value pairs
            gene_id = None
            gene_name = None
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
            if gene_id:
                gene_map[gene_id] = gene_name if gene_name else gene_id
    return gene_map


def extract_gene_ids_from_isoform_str(isoform_str: str, gene_map: Dict[str, str]) -> Tuple[str, str]:
    """
    Extract gene IDs and names from isoform_ids string.
    
    Args:
        isoform_str: comma-separated string of isoform_id_gene_id pairs
        gene_map: dict mapping gene_id -> gene_name
    
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
            gene_id = parts[-1]  # gene_id is after last underscore
            if gene_id.startswith('ENS'):
                ens_gene_ids.append(gene_id)
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


def annotate_isoform_events(input_df: pd.DataFrame, gene_map: Dict[str, str]) -> pd.DataFrame:
    """
    Annotate isoform events dataframe with gene IDs and names.
    
    Args:
        input_df: pandas DataFrame with 'isoform_ids' column
        gene_map: dict mapping gene_id -> gene_name
    
    Returns:
        pandas DataFrame with added 'gene_id' and 'gene_name' columns
    """
    if 'isoform_ids' not in input_df.columns:
        raise ValueError("Input DataFrame must contain 'isoform_ids' column")
    
    gene_ids = []
    gene_names = []
    warning_count = 0
    
    for idx, isoform_str in enumerate(input_df['isoform_ids']):
        gene_id_str, gene_name_str = extract_gene_ids_from_isoform_str(isoform_str, gene_map)
        
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


def annotate_single_file(input_file: str, gene_map: dict, output_file: str = "", progress=None) -> str:
    """
    Annotate a single TSV file and return output path.
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
    df = annotate_isoform_events(df, gene_map)
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
    try:
        import matplotlib.pyplot as plt
        import numpy as np
    except ImportError:
        import sys
        print("Warning: matplotlib not installed. Skipping plot generation.", file=sys.stderr)
        return
    
    # Sort types in order: ES, Alt5, Alt3, IR (if present)
    type_order = ['ES', 'Alt5', 'Alt3', 'IR']
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
    ax.barh(y, inclusion_counts, height, label='Inclusion', color='#1f77b4')
    ax.barh(y, exclusion_counts, height, left=inclusion_counts, label='Exclusion', color='#ff7f0e')
    
    ax.set_xlabel('Count')
    ax.set_ylabel('Splicing Type')
    ax.set_title('Inclusion/Exclusion Events by Splicing Type')
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
    
    plot_path = f"{output_prefix}_splicing_counts.png"
    # Ensure output directory exists
    from pathlib import Path
    plot_path_obj = Path(plot_path)
    plot_path_obj.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(plot_path, dpi=300)
    plt.close()
    print(f"Plot saved to: {plot_path}")