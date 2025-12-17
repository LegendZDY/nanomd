import typer
import subprocess
from pathlib import Path
from rich.progress import track
from typing_extensions import Annotated
from basebio import minimap2

app = typer.Typer()

@app.command()
def isoformAS(
    input: Annotated[str, typer.Option("--input", "-i", help="Input fastq files.")],
    reference: Annotated[str, typer.Option("--reference", "-r", help="reference genome path.")],
    prefix: Annotated[str, typer.Option("--prefix", "-p", help="Prefix for output files.")],
    tool: Annotated[str, typer.Option("--tool", help="minimap2.")]="minimap2",
    parms: Annotated[str, typer.Option("--parms", help="minimap2 parameters for mapping.")]="--secondary=no --cs -a",
    threads: Annotated[int, typer.Option(help="Number of threads.")]=4,
    ):
    """
    Alternative splicing analysis.
    """
    minimap2(input, reference, prefix, tool, parms, threads)