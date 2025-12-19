import time
import typer
from typing_extensions import Annotated
from rich.progress import Progress, SpinnerColumn, TextColumn
from basebio import minimap2
from basebio import check_path_exists
from ..utils.quantify import salmon_map, salmon_quantify


app = typer.Typer()

@app.command()
def count(
    input: Annotated[str, typer.Option("--input", "-i", help="Input fastq files.")],
    reference: Annotated[str, typer.Option("--reference", "-r", help="reference transcripts path.")],
    output: Annotated[str, typer.Option("--output", "-o", help="Output file path.")]=".",
    prefix: Annotated[str, typer.Option("--prefix", "-p", help="Prefix for output files.")]="prefix",
    ):
    """
    Mapping of nanopore reads to transcripts reference.
    """
    
    with Progress(
        SpinnerColumn(),
        TextColumn("[progress.description]{task.description}"),
        transient=True,
    ) as progress:
        try:
            progress.add_task(description="map transcripts reference...", total=None)
            start=time.time()
            
            progress.add_task(description="Mapping reads to reference...", total=None)
            output_bam=f"{output}/{prefix}_transcripts.bam"
            if not check_path_exists(output_bam):
                salmon_map(input, reference, output_bam)
            progress.add_task(description="Mapping reads to reference Done", total=None)

            output_quant=f"{output}/{prefix}_quant"
            progress.add_task(description="Quantifying reads ...", total=None)
            if not check_path_exists(output_quant):
                salmon_quantify(output_bam, reference, output_quant)
            progress.add_task(description="Quantifying reads Done", total=None)
            end=time.time()

            time_cost=f"{(end - start) // 3600}h{((end - start) % 3600) // 60}m{(end - start) % 60:.2f}s"
            print(f"map transcripts reference Done, time cost: {time_cost}")
            progress.add_task(description=f"map transcripts reference Done, time cost: {time_cost}", total=None)
        except Exception as e:
            print(f"Error: {e}")
            progress.add_task(description="map transcripts reference Failed", total=None)
            exit(1)
    