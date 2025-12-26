import time, os
import typer
from typing_extensions import Annotated
from rich.progress import Progress, SpinnerColumn, TextColumn
from pathlib import Path
from basebio import check_path_exists, run_command
from ..utils.modtools import split_mod
from ..utils.abs_position import gene_feature_distance_calculator
from ..utils.modifications import form_reads_get_modifications

app = typer.Typer()

@app.command()
def detectMod(
    input: Annotated[str, typer.Option("--input", "-i", help="Input fastq files.")],
    sam: Annotated[str, typer.Option("--sam", "-s", help="mapping sam/bam file.")],
    bed: Annotated[str, typer.Option("--bed", "-b", help="bed file for transcripts sites.")],
    regions: Annotated[str, typer.Option("--regions", "-r", help="regions file for gene feature distance calculation.")],
    output: Annotated[str, typer.Option("--output", "-o", help="Output file path.")]=".",
    prefix: Annotated[str, typer.Option("--prefix", "-p", help="Prefix for output files.")]="prefix",
    pvalue: Annotated[float, typer.Option("--pvalue", help="pvalue cutoff for modification sites.")]=0.98,
    ):
    """
    Detect modification sites in input fastq files.
    """
    
    with Progress(
        SpinnerColumn(),
        TextColumn("[progress.description]{task.description}"),
        transient=True,
    ) as progress:
        progress.add_task(description="Detecting modification start...", total=None)
        start=time.time()
        plot_script = Path(__file__).parent.parent / "scripts"
        WKD = os.getcwd()
        
        output_file=f"{output}/{prefix}.bed"
        progress.add_task(description="Getting modification from fq files...", total=None)
        if not check_path_exists(output_file):
            mod = form_reads_get_modifications(input, sam, bed, output_file, pvalue)
            mod.get_mod_position_with_sam()
        progress.add_task(description=f"Getting modification from fq files Done", total=None)

        progress.add_task(description="Splitting modification sites...", total=None)
        if not check_path_exists(f"{output}/{prefix}_m6A.bed"):
            split_mod(output_file, prefix)
        progress.add_task(description="Splitting modification sites Done", total=None)

        progress.add_task(description="Calculating absolute distance...", total=None)
        bed_file_list = [f"{output}/{prefix}_m6A.bed", f"{output}/{prefix}_m5C.bed", 
                         f"{output}/{prefix}_psi.bed", f"{output}/{prefix}_AtoI.bed"]
        for bed_file in bed_file_list:
            gfc_output = bed_file.replace(".bed", "_abs_dist.txt")
            if not check_path_exists(gfc_output):
                gfc = gene_feature_distance_calculator(bed_file, regions, gfc_output)
                gfc.process_bed_file()
        progress.add_task(description="Calculating absolute distance Done", total=None)
        
        end=time.time()
        time_cost=f"{(end - start) // 3600}h{((end - start) % 3600) // 60}m{(end - start) % 60:.2f}s"
        print(f"Detecting modification sites Done, time cost: {time_cost}")
        progress.add_task(description=f"Detecting modification sites Done, time cost: {time_cost}", total=None)