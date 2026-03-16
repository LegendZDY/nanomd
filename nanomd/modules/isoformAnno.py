import time
import typer
import pandas as pd
from typing_extensions import Annotated
from rich.progress import Progress, SpinnerColumn, TextColumn
from pathlib import Path
import glob
import sys
import os
from ..utils.annotation_tools import parse_gtf, annotate_isoform_events, annotate_single_file, count_inclusion_exclusion, plot_splicing_counts

app = typer.Typer()

@app.command()
def isoformAnno(
    input: Annotated[str, typer.Option("--input", "-i", help="Input TSV file or glob pattern (e.g., diffsplice*.tsv) of alternative splicing events.")],
    gtf: Annotated[str, typer.Option("--gtf", "-g", help="GTF annotation file.")],
    output: Annotated[str, typer.Option("--output", "-o", help="Output annotated TSV file (if single input) or output prefix (if multiple files).")]="",
    plot: Annotated[bool, typer.Option("--plot", help="Generate stacked bar chart of inclusion/exclusion counts.")] = False,
    ):
    """
    Annotate alternative splicing events with gene ID and gene name.
    Supports single file or glob pattern for multiple files.
    """
    start_time = time.time()
    
    # Determine if input is a glob pattern
    input_files = glob.glob(input)
    if not input_files:
        print(f"Error: No files found matching pattern: {input}", file=sys.stderr)
        exit(1)
    
    # Filter to only TSV files that match diffsplice pattern
    splicing_types = ['es', 'alt5', 'alt3', 'ir']
    filtered_files = []
    for f in input_files:
        fname = os.path.basename(f).lower()
        if fname.endswith('.tsv') and any(st in fname for st in splicing_types):
            filtered_files.append(f)
    
    if not filtered_files:
        print(f"Error: No diffsplice TSV files found in: {input_files}", file=sys.stderr)
        exit(1)
    
    # Sort files by splicing type for consistent ordering
    filtered_files.sort(key=lambda x: os.path.basename(x).lower())
    
    with Progress(
        SpinnerColumn(),
        TextColumn("[progress.description]{task.description}"),
        transient=True,
    ) as progress:
        try:
            progress.add_task(description="Loading GTF annotation...", total=None)
            gene_map = parse_gtf(gtf)
            progress.add_task(description=f"Loaded {len(gene_map)} gene mappings", total=None)
            
            # Process each file
            output_files = []
            counts_dict = {}  # {splicing_type: (inclusion, exclusion)}
            
            for input_file in filtered_files:
                fname = os.path.basename(input_file)
                # Determine splicing type from filename
                splicing_type = 'unknown'
                if 'es' in fname.lower():
                    splicing_type = 'ES'
                elif 'alt5' in fname.lower():
                    splicing_type = 'Alt5'
                elif 'alt3' in fname.lower():
                    splicing_type = 'Alt3'
                elif 'ir' in fname.lower():
                    splicing_type = 'IR'
                
                # Determine output file name
                if len(filtered_files) == 1 and output:
                    # Single file mode with explicit output name
                    output_file = output
                else:
                    # Multiple files or default output name
                    input_path = Path(input_file)
                    if output:
                        # Use output as prefix
                        output_file = f"{output}_{splicing_type}.annotated.tsv"
                    else:
                        output_file = str(input_path.with_suffix('.annotated.tsv'))
                
                # Annotate the file
                annotated_path = annotate_single_file(input_file, gene_map, output_file, progress)
                output_files.append(annotated_path)
                
                # Count inclusion/exclusion
                df = pd.read_csv(input_file, sep='\t')
                inc, exc = count_inclusion_exclusion(df)
                counts_dict[splicing_type] = (inc, exc)
                
                progress.add_task(description=f"{splicing_type}: {inc} inclusion, {exc} exclusion events", total=None)
            
            # Generate plot if requested
            if plot:
                progress.add_task(description="Generating plot...", total=None)
                if output:
                    plot_prefix = output
                else:
                    # Use first input file's directory and base name
                    first_path = Path(filtered_files[0])
                    plot_prefix = str(first_path.parent / first_path.stem)
                plot_splicing_counts(counts_dict, plot_prefix)
            
            end_time = time.time()
            time_cost = f"{(end_time - start_time) // 3600}h{((end_time - start_time) % 3600) // 60}m{(end_time - start_time) % 60:.2f}s"
            
            print(f"\nAnnotation completed, time cost: {time_cost}")
            print(f"Processed {len(filtered_files)} files:")
            for splicing_type, (inc, exc) in counts_dict.items():
                print(f"  {splicing_type}: {inc} inclusion, {exc} exclusion events")
            print(f"Output files:")
            for out_file in output_files:
                print(f"  {out_file}")
            
        except Exception as e:
            print(f"Error: {e}", file=sys.stderr)
            progress.add_task(description="Annotation failed", total=None)
            exit(1)