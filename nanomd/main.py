import typer
from .modules.gene import gene
from .modules.count import count
from .modules.matrix import matrix
from .modules.polyA import polyA
from .modules.detectMod import detectMod
from .modules.isoformAS import isoformAS
from .modules.nascentRNA import nascentRNA

app = typer.Typer(add_completion=False)

@app.callback()
def callback():
    """
    nanoMD(Nanopore direct RNA sequencing Multi-dimensional analysis) was developed to synchronously analyze the changes in m6A, m5C, psi, AtoI modification sites, genes, isoforms, alternative splicing events, nascent RNA and polyA tail in direct RNA sequencing data.
    """

app.command(name="gene")(gene)
app.command(name="count")(count)
app.command(name="polyA")(polyA)
app.command(name="matrix")(matrix)
app.command(name="isoformAS")(isoformAS)
app.command(name="detectMod")(detectMod)
app.command(name="nascentRNA")(nascentRNA)


if __name__ == "__main__":
    app()