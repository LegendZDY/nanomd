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