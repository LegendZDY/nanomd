import time
import typer
import pandas as pd
from typing_extensions import Annotated
from rich.progress import Progress, SpinnerColumn, TextColumn
from pathlib import Path
from ..utils.annotation_tools import parse_gtf, annotate_isoform_events

app = typer.Typer()

@app.command()
def isoformAnno(
    input: Annotated[str, typer.Option("--input", "-i", help="Input TSV file of alternative splicing events.")],
    gtf: Annotated[str, typer.Option("--gtf", "-g", help="GTF annotation file.")],
    output: Annotated[str, typer.Option("--output", "-o", help="Output annotated TSV file.")]="",
    ):
    """
    Annotate alternative splicing events with gene ID and gene name.
    """
    if not output:
        input_path = Path(input)
        output = str(input_path.with_suffix('.annotated.tsv'))
    
    with Progress(
        SpinnerColumn(),
        TextColumn("[progress.description]{task.description}"),
        transient=True,
    ) as progress:
        try:
            progress.add_task(description="Loading GTF annotation...", total=None)
            start = time.time()
            gene_map = parse_gtf(gtf)
            progress.add_task(description=f"Loaded {len(gene_map)} gene mappings", total=None)
            
            progress.add_task(description="Reading input TSV...", total=None)
            df = pd.read_csv(input, sep='\t')
            
            progress.add_task(description="Annotating isoform events...", total=None)
            df = annotate_isoform_events(df, gene_map)
            
            progress.add_task(description="Writing output...", total=None)
            df.to_csv(output, sep='\t', index=False)
            
            end = time.time()
            time_cost = f"{(end - start) // 3600}h{((end - start) % 3600) // 60}m{(end - start) % 60:.2f}s"
            print(f"Annotation completed, time cost: {time_cost}")
            print(f"Output saved to: {output}")
            progress.add_task(description=f"Annotation completed, time cost: {time_cost}", total=None)
        except Exception as e:
            print(f"Error: {e}")
            progress.add_task(description="Annotation failed", total=None)
            exit(1)