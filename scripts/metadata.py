import tempfile
import tarfile
import subprocess
from urllib.request import urlretrieve
import ncbitax2lin  # just to make sure it's there
import pandas as pd
import os


def main():
    os.makedirs("intermediate_files", exist_ok=True)
    metadata = download_metadata()
    metadata.to_csv("intermediate_files/metadata.csv", index=False)


def download_metadata():
    with tempfile.TemporaryDirectory(prefix="cocopye_", dir="intermediate_files") as tmpdir:
        print("  - Downloading accession metadata")
        urlretrieve(
            "https://ftp.ncbi.nlm.nih.gov/genomes/refseq/assembly_summary_refseq.txt",
            os.path.join(tmpdir, "accession_metadata.tsv")
        )

        print("  - Downloading taxdump")
        urlretrieve(
            "https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz",
            os.path.join(tmpdir, "taxdump.tar.gz")
        )

        print("  - Extracting taxdump")
        file = tarfile.open(os.path.join(tmpdir, "taxdump.tar.gz"))
        file.extractall(os.path.join(tmpdir, "taxdump"))
        file.close()

        print("  - Running ncbitax2lin")
        tax2lin_process = subprocess.Popen(
            [
                "ncbitax2lin",
                "--nodes-file", "taxdump/nodes.dmp",
                "--names-file", "taxdump/names.dmp",
                "--output", "tax2lin.csv.gz"
            ],
            cwd=tmpdir
        )
        tax2lin_process.wait()

        print("  - Reading data")
        acc2tax = pd.read_csv(os.path.join(tmpdir, "accession_metadata.tsv"), sep="\t", skiprows=1)
        tax2lin = pd.read_csv(os.path.join(tmpdir, "tax2lin.csv.gz"), sep=",", compression='gzip')

        print("  - Merging data")
        merged = acc2tax.merge(tax2lin, left_on="taxid", right_on="tax_id")
        merged = merged[["#assembly_accession", "superkingdom", "phylum", "class", "order", "family", "genus", "species"]]
        merged = merged.rename(columns={"#assembly_accession": "sequence"})

        return merged 


if __name__ == "__main__":
    main()
