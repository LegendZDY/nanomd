import time, os
import typer
from typing_extensions import Annotated
from rich.progress import Progress, SpinnerColumn, TextColumn
from pathlib import Path
from basebio import check_path_exists, run_command
from ..utils.quantify import matrix_generate
from ..utils.ployA_tools import ployA_matrix_generate

app = typer.Typer()

@app.command()
def matrix(
    inputs: Annotated[str, typer.Option("--inputs", "-i", help="Regular matching pattern for salmon results directories or ployA files, such as 'path/to/*_quant' or 'path/to/*ployA.tsv'.")],
    control_names: Annotated[str, typer.Option("--control_names", "-c", help="Control sample dir name or ployA file name separated by comma, such as 'control1_quant,control2_quant,control3_quant' or 'control1_ployA.tsv,control2_ployA.tsv,control3_ployA.tsv'.")],
    prefix: Annotated[str, typer.Option("--prefix", "-p", help="Prefix for output file.")] = "prefix",
    species: Annotated[str, typer.Option("--species", "-s", help="Species name. [human|mouse|rat]")] = "human",
    type: Annotated[str, typer.Option("--type", "-t", help="Input count file type, either 'salmon', 'ployA', 'meta'.")] = "salmon",
    docker: Annotated[bool, typer.Option("--docker", help="Whether to run in docker container to plots.")] = False,
    ):
    """
    Generate matrix with salmon results or ployA files.
    """
    
    with Progress(
        SpinnerColumn(),
        TextColumn("[progress.description]{task.description}"),
        transient=True,
    ) as progress:
        try:
            start=time.time()

            task1 = progress.add_task(description="Generate count matrix and plot matrix...", total=100)
            plot_script = Path(__file__).parent.parent / "scripts"
            split_num = len(control_names.split(","))
            WKD = os.getcwd()

            if type == "salmon":
                output_count = "matrix_count.tsv"
                if not check_path_exists(output_count):
                    matrix_generate(input_dirs=inputs, control_dirs_names=control_names, output=output_count,
                                    count_type="NumReads")
                output_tpm = "matrix_tpm.tsv"
                if not check_path_exists(output_tpm):
                    matrix_generate(input_dirs=inputs, control_dirs_names=control_names, output=output_tpm,
                                    count_type="TPM")
                output_plot = f"matrix_plots"
                if not check_path_exists(output_plot):
                    if docker:
                        run_command(f"docker run -v {WKD}:/output -v {plot_script}:/scripts -w /output legendzdy/rbase:1.0.0 Rscript /scripts/isoform.R -i /output/{output_count} -o /output/{output_plot} -p {prefix} -s {species} -t {type} -n {split_num}".split())
                    else:
                        run_command(f"Rscript {plot_script}/isoform.R -i {WKD}/{output_count} -o {WKD}/{output_plot} -p {prefix} -s {species} -t {type} -n {split_num}")
            elif type == "ployA":
                output_count = "matrix_ployA.tsv"
                output_lengths = "ployA_lengths.tsv"
                if not check_path_exists(output_count) or not check_path_exists(output_lengths):
                    ployA_matrix_generate(input_files=inputs, control_filess_names=control_names, output_matrix=output_count, output_lengths=output_lengths)
                output_plot = f"matrix_plots"
                if not check_path_exists(output_plot):
                    if docker:
                        run_command(f"docker run -v {WKD}:/output -v {plot_script}:/scripts -w /output legendzdy/rbase:1.0.0 Rscript /scripts/isoform.R -i /output/{output_count} -o /output/{output_plot} -p {prefix} -s {species} -t {type} -n {split_num}".split())
                    else:
                        run_command(f"Rscript {plot_script}/isoform.R -i {WKD}/{output_count} -o {WKD}/{output_plot} -p {prefix} -s {species} -t {type} -n {split_num}")
            elif type == "meta":
                output_plot = f"meta_plots"
                for fileType in ["m6A", "m5C", "psi", "AtoI"]:
                    output_plot = f"{prefix}_{fileType}_metagene.pdf"
                    if not check_path_exists(output_plot):
                        if docker:
                            run_command(f"docker run -v {WKD}:/output -v {plot_script}:/scripts -w /output legendzdy/rbase:1.0.0 Rscript /scripts/metaplot.R -i /output -o /output/{output_plot} -p {prefix} -t {fileType}".split())
                        else:
                            run_command(f"Rscript {plot_script}/metaplot.R -i {WKD} -o {WKD}/{output_plot} -p {prefix} -t {fileType}")
            else:
                print(f"Error: Unknown count type {type}")
                exit(1)
            progress.update(task1, advance=100)
            end=time.time()
            time_cost=f"{(end - start) // 3600}h{((end - start) % 3600) // 60}m{(end - start) % 60:.2f}s"
            print(f"Generate count matrix and plot matrix Done, time cost: {time_cost}")
            progress.console.print(f"Generate count matrix and plot matrix Done, time cost: {time_cost}")
        except Exception as e:
            print(f"Error: {e}")
            progress.console.print("Generate count matrix and plot matrix Failed")
            exit(1)