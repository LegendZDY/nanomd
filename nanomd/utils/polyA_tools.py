# -*- coding: utf-8 -*-
__author__ = "legendzdy@dingtalk.com"
"""
Author: legendzdy@dingtalk.com
Data: 20250612
Description: This file contains the functions for converting pod5 files to fast5 files and indexing the fastq file using nanopolish index.
"""
from concurrent.futures import Future, ProcessPoolExecutor, as_completed
from pathlib import Path
from typing import Dict, List
from more_itertools import chunked

import pysam, os, glob
import pod5 as p5
import pandas as pd
from pod5.tools.utils import DEFAULT_THREADS, collect_inputs, limit_threads
from pod5.tools.pod5_convert_to_fast5 import * # type: ignore
from basebio import run_command

def convert_to_fast5_with_summary_file(
    inputs: List[Path],
    output: Path,
    summary_file: str,
    fastq_name: str,
    recursive: bool = False,
    threads: int = DEFAULT_THREADS,
    force_overwrite: bool = False,
    file_read_count: int = 4000,
):
    if output.exists() and not output.is_dir():
        raise FileExistsError("Cannot output to a file")

    threads = limit_threads(threads)
    columns = ["filename_fastq", "filename_fast5", "read_id", "run_id", "channel", "mux", "start_time", "duration", "num_events", "passes_filtering", "template_start", "num_events_template", "template_duration", "sequence_length_template", "mean_qscore_template", "strand_score_template", "median_template", "mad_template", "pore_type", "experiment_id", "sample_id", "end_reason", "alias", "type", "barcode_arrangement", "barcode_full_arrangement", "barcode_kit", "barcode_variant", "barcode_score", "barcode_front_id", "barcode_front_score", "barcode_front_refseq", "barcode_front_foundseq", "barcode_front_foundseq_length", "barcode_front_begin_index", "barcode_rear_id", "barcode_rear_score", "barcode_rear_refseq", "barcode_rear_foundseq", "barcode_rear_foundseq_length", "barcode_rear_end_index"]

    # 数据行（单行列表）
    data_row = ["test", "885", "1", "56.132750", "0.756500", "605", "TRUE", "56.154000", "588", "0.735250", "153", "8.314826", "0.000000", "94.303551", "8.406906", "not_set", "20220426-CDS0204-PAK00515", "no_sample", "signal_positive", "barcode11", "test_sample", "barcode11", "LWB11_var1", "LWB", "var1", "92.750000", "BC11_FWD", "89.250000", "CCGTGACGTT", "51", "34", "BC11_REV", "92.750000", "GCAATAA", "GCACATC", "52", "20"]
    
    with open(summary_file, "w") as fo:
        fo.write("\t".join(columns) + "\n")
        with ProcessPoolExecutor(max_workers=threads) as executor:
            total_reads = 0
            futures: Dict[Future, Path] = {}
            # Enumerate over input pod5 files
            for input_idx, source in enumerate(
                collect_inputs(inputs, recursive, "*.pod5", threads=threads)
            ):
                # Open the inputs to read the read ids
                with p5.Reader(source) as reader:
                    for chunk_idx, read_ids in enumerate(
                        chunked(reader.read_ids, file_read_count)
                    ):
                        dest = (
                            output / f"{source.stem}.{chunk_idx}_{input_idx}.fast5"
                        ).resolve()

                        if dest.exists() and not force_overwrite:
                            raise FileExistsError(
                                "Output path points to an existing file and --force-overwrite not set"
                            )

                        kwargs = {
                            "source": source,
                            "dest": dest,
                            "read_ids": read_ids,
                        }
                        futures[executor.submit(convert_pod5_to_fast5, **kwargs)] = dest  # type: ignore

                        for read_id in read_ids:
                            raw = [fastq_name, f"{source.stem}.{chunk_idx}_{input_idx}.fast5", read_id] + data_row
                            fo.write("\t".join(raw) + "\n")

                    total_reads += len(reader.read_ids)

        print(f"Converting pod5s into {len(futures)} fast5 files. Please wait...")

        status = StatusMonitor(file_count=len(inputs))
        status.increment(files_started=len(inputs), read_count=total_reads)

        for idx, future in enumerate(as_completed(futures)):
            (reads_converted, samples_converted) = future.result()

            status.increment(
                files_ended=1,
                sample_count=samples_converted,
                reads_processed=reads_converted,
            )
            status.print_status()

        status.print_status(force=True)

        print("Conversion complete")

def index_fastq(fast5_dir, summary_file, fastq):
    """
    Index the fastq file using nanopolish index.
    
    Args:
        fast5_dir (str): The directory where the fast5 files are located.
        summary_file (str): The path to the sequencing summary file.
        fastq (str): The path to the fastq file.
    
    Examples:
        index_fastq("fast5_dir", "summary.txt", "fastq.gz")
    """

    cmd = f"nanopolish index --directory={fast5_dir} --sequencing-summary={summary_file} {fastq}".split()
    run_command(cmd)

def detect_polyA(fastq, bam, transcriptome, output_ploya, threads=8):
    """
    Detect polyA using nanopolish polyA.

    Args:
        fastq (str): The path to the fastq file.
        bam (str): The path to the bam file.
        transcriptome (str): The path to the transcriptome file.
        output_ploya (str): The path to the output polyA file.
        threads (int): The number of threads to use.

    Examples:
        detect_polyA("fastq.gz", "bam", "transcriptome.fa", "polyA.tsv", threads=8)
    """

    cmd = f"nanopolish polyA --threads={threads} --reads={fastq} --bam={bam} --genome={transcriptome} > {output_ploya}"
    run_command(cmd, use_shell=True)

class polyADetector:
    """
    polyA detector.

    Args:
        bam_path (str): The path to the bam file.
        output_path (str): The path to the output file.
        min_a_length (int): The minimum length of polyA.
        max_non_a (int): The maximum number of non-A bases allowed in polyA.
    
    Examples:
        detector = polyADetector("bam", "polyA.tsv", min_a_length=12, max_non_a=3)
        detector.analyze()
    """
    def __init__(self, bam_path, output_path, min_a_length=12, max_non_a=3):
        self.bam_path = bam_path
        self.output_path = output_path
        self.min_a_length = min_a_length
        self.max_non_a = max_non_a
    
    @staticmethod
    def reverse_complement(seq):
        """
        Get the reverse complement of a sequence.
        """
        complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'U': 'A',
                    'a': 't', 't': 'a', 'g': 'c', 'c': 'g', 'u': 'a',
                    'N': 'N', 'n': 'n'}
        return ''.join(complement.get(base, base) for base in reversed(seq)) # type: ignore
    
    def find_longest_polyA(self, seq):
        """
        Use sliding window algorithm to find the longest consecutive A sequence, allowing a few non-A bases.
        """
        if not seq:
            return "", 0, 0, 0.0
        
        n = len(seq)
        left = 0
        non_a_count = 0
        max_length = 0
        a_count_in_max = 0
        best_start = 0
        best_end = -1
        
        for right in range(n):
            if seq[right].upper() != 'A':
                non_a_count += 1
            
            while non_a_count > self.max_non_a:
                if seq[left].upper() != 'A':
                    non_a_count -= 1
                left += 1
            
            current_length = right - left + 1
            current_a_count = current_length - non_a_count

            if current_length > max_length:
                max_length = current_length
                a_count_in_max = current_a_count
                best_start = left
                best_end = right
        
        a_ratio = a_count_in_max / max_length if max_length > 0 else 0.0
        
        polyA_seq = seq[best_start:best_end+1] if max_length > 0 else ""
        
        return polyA_seq, max_length, a_count_in_max, a_ratio
    
    def process_read(self, read):
        if read.is_unmapped or not read.cigartuples:
            return None
        
        read_name = read.query_name
        seq = read.query_sequence
        if not seq:
            return None
            
        flag = read.flag
        cigar = read.cigartuples
        
        if flag == 0:
            target_op = cigar[-1]
        elif flag == 16:
            target_op = cigar[0]
        else:
            return None
        
        op_type, op_len = target_op
        
        if op_type != 4:
            return None
        
        if flag == 0:
            target_seq = seq[-op_len:] if op_len <= len(seq) else seq
        elif flag == 16:
            target_seq = self.reverse_complement(seq[:op_len]) if op_len <= len(seq) else seq
        else:
            return None
        
        polyA_seq, total_length, a_count, a_ratio = self.find_longest_polyA(target_seq)
        
        has_polyA = a_count >= self.min_a_length and a_ratio >= 0.75
        
        return {
            "read_name": read_name,
            "ref_name": read.reference_name,
            "strand": "+" if flag == 0 else "-",
            "polyA_seq": polyA_seq,
            "polyA_region_length": total_length,
            "a_count": a_count,
            "a_ratio": a_ratio,
            "has_polyA": "Yes" if has_polyA else "No"
        }
    
    def analyze(self):
        with open(self.output_path, 'w') as fout:
            fout.write("readName\trefName\tstrand\tpolyASeq\tpolyALength\tACount\tARatio\tHaspolyA\n")
            with pysam.AlignmentFile(self.bam_path, "rb") as bam:
                for read in bam:
                    result = self.process_read(read)
                    if result:
                        if result["has_polyA"] == "Yes":
                            fout.write("\t".join([
                                result["read_name"],
                                result["ref_name"],
                                result["strand"],
                                result["polyA_seq"],
                                str(result["polyA_region_length"]),
                                str(result["a_count"]),
                                f"{result['a_ratio']:.3f}",
                                result["has_polyA"]
                            ]) + "\n")

def drop_outliers(df, dep_col, handle_nan='keep'):
    """
    Use IQR method to drop outliers from DataFrame

    Args:
        df: pandas DataFrame
        dep_col: str, the column name for detecting outliers
        handle_nan: str, the method for handling NaN values ('keep','remove', 'fill')
    
    Returns:
        The filtered DataFrame
    """

    # Validate column exists
    if dep_col not in df.columns:
        raise ValueError(f"Column '{dep_col}' not found in DataFrame")

    # Handle NaN values
    if handle_nan == 'remove':
        df = df.dropna(subset=[dep_col])
    elif handle_nan == 'fill':
        median_val = df[dep_col].median()
        df[dep_col] = df[dep_col].fillna(median_val)

    # Default 'keep' does not do anything
    # Calculate quartiles and IQR
    a = df[dep_col]
    q1 = a.quantile(0.25)
    q3 = a.quantile(0.75)
    iqr = q3 - q1

    # Calculate normal value range
    bottom = q1 - 1.5 * iqr
    top = q3 + 1.5 * iqr

    # Create boolean mask
    mask = (a >= bottom) & (a <= top)
    
    return df[mask]

def polyA_matrix_generate(input_files, control_filess_names, output_matrix, output_lengths):
    """
    Generate polyA matrix and polyA lengths file.
    """
    if isinstance(input_files, str):
        sorted_files = sorted(glob.glob(input_files))
        if not sorted_files:
            raise ValueError(f"No matching file found for {input_files}")
        input_path = os.path.dirname(sorted_files[0])
        control_files = [os.path.join(input_path, control_name) for control_name in control_filess_names.split(',')]
        for control_file in control_files:
            if control_file in sorted_files:
                sorted_files.remove(control_file)   
            else:
                raise ValueError(f"Control file {control_file} not found in {input_files}")
        input_files = control_files + sorted_files
    elif not input_files:
        raise ValueError("input_files should not be empty")
    
    # Initialize a dictionary to store all counts
    all_counts = {}
    sample_names = []
    # 创建一个空的DataFrame存储所有结果
    all_polyA_lengths = pd.DataFrame()
    
    # Process each sample file
    for sample_file in input_files:
        sample_name = os.path.basename(sample_file).split('.')[0].replace('_polyA', '')
        sample_names.append(sample_name)
        condition = sample_name
        
        # Read the TSV file
        df = pd.read_csv(sample_file, sep='\t')
        
        # Count occurrences of each refName
        ref_counts = df['refName'].value_counts().to_dict()
        all_counts[sample_name] = ref_counts
        # polyA lengths
        df_temp = df[['polyALength']].copy()
        df_temp['condition'] = condition
        df_temp = df_temp.rename(columns={'polyALength': 'measurement'})
        all_polyA_lengths = pd.concat([all_polyA_lengths, df_temp], ignore_index=True)
    
    # Create a set of all unique refNames
    all_refnames = set()
    for counts in all_counts.values():
        all_refnames.update(counts.keys())
    
    # Create the matrix DataFrame
    matrix = pd.DataFrame(index=list(all_refnames), columns=sample_names)
    matrix = matrix.fillna(0)  # Fill with zeros initially
    # Populate the matrix with counts
    for sample, counts in all_counts.items():
        for refname, count in counts.items():
            matrix.at[refname, sample] = count
    # Save the matrix to output file
    matrix.to_csv(output_matrix, sep='\t')

    # polyA lengths
    all_polyA_lengths = all_polyA_lengths[['condition', 'measurement']]
    all_polyA_lengths = drop_outliers(all_polyA_lengths, 'measurement')
    all_polyA_lengths.to_csv(output_lengths, sep='\t', index=False)