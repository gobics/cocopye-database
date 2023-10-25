import os
import datetime
import gzip
import shutil
import tempfile
import tarfile
import subprocess
from urllib.request import urlretrieve

import pandas as pd
import numpy as np
import ncbi_genome_download as ngd
from cocopye.pfam import count_pfams
from cluster import cluster


# === Change this according to your system and requirements. ===============

PFAM_VERSION = ["28", "24"]  # Can also contain both version 24 and 28
MIN_COMP = 0.95
MAX_CONT = 0.05

UPROC_ORF = "/home/niklas/.local/share/cocopye/uproc/bin/uproc-orf"
UPROC_PROT = "/home/niklas/.local/share/cocopye/uproc/bin/uproc-prot"
PFAM_DIR = "/home/niklas/.local/share/cocopye/pfam_db"
MODEL_DIR = "/home/niklas/.local/share/cocopye/model"

# ==========================================================================


def main():
    # Download sequences
    #log("Downloading sequences. This may take a while (up to several days)")
    #ngd.download(assembly_levels="complete,chromosome", file_formats="fasta", output="intermediate_files", groups="archaea,bacteria")

    # Extract sequences
    #log("Extracting archaea files")
    #extract_files("intermediate_files/refseq/archaea", "intermediate_files/fasta_files")
    #log("Extracting bacteria files")
    #extract_files("intermediate_files/refseq/bacteria", "intermediate_files/fasta_files")

    # Load metadata
    log("Loading metadata")
    metadata = pd.read_csv("intermediate_files/metadata.csv")

    archaea = list(metadata.loc[metadata["superkingdom"] == "Archaea"]["sequence"])
    bacteria = list(metadata.loc[metadata["superkingdom"] == "Bacteria"]["sequence"])

    for version in PFAM_VERSION:
        # Count Pfams
        log("Generating Pfam count matrix (Pfam v" + version + ")")
        counts, sequences = count_pfams(UPROC_ORF, UPROC_PROT, os.path.join(PFAM_DIR, version), MODEL_DIR, "intermediate_files/fasta_files")[:2]
        counts = counts[:,1:]
        pfam_counts = pd.DataFrame(counts, index=sequences, columns=[f"PF{'0' * (5-len(str(i)))}{i}" for i in range(1, 17127)])
        os.makedirs(os.path.join("intermediate_files", version), exist_ok=True)
        pfam_counts.to_csv(os.path.join("intermediate_files", version, "pfam_counts.csv"))

        # Cluster sequences
        log("Clustering sequences")
        reps = cluster(pfam_counts)

        # Calculate universal markers
        log("Calculating universal markers based on clustered sequences")
        universal_bac = universal_markers(reps, bacteria)
        universal_arc = universal_markers(reps, archaea)

        # Filter sequences
        log("Filtering all sequences based on universal markers")
        filtered_bac = filter_sequences(pfam_counts, bacteria, universal_bac)
        filtered_arc = filter_sequences(pfam_counts, archaea, universal_arc)
        filtered_seq = pd.concat([filtered_bac, filtered_arc])

        # Cluster again
        log("Clustering filtered sequences")
        final_sequences = cluster(filtered_seq)
        final_sequences.to_csv(os.path.join("intermediate_files", version, "filtered.csv"))

        # Move everything to final output folder
        log("Move everything to final output folder")
        final_seq_ids = list(final_sequences.index)
        final_metadata = metadata.loc[metadata["sequence"].isin(final_seq_ids)]
        os.makedirs(os.path.join("output", version, "fasta"), exist_ok=True)
        final_metadata.to_csv(os.path.join("output", version, "metadata.csv"))
        for seq_id in final_seq_ids:
            shutil.copy(
                os.path.join("intermediate_files", "fasta_files", seq_id + ".fna"),
                os.path.join("output", version, "fasta", seq_id + ".fna")
            )


def log(msg):
    print("[" + str(datetime.datetime.now()) + "] " + msg)


def extract_files(input_path, output_path):
    os.makedirs(output_path, exist_ok=True)

    dirs = [p for p in os.listdir(input_path) if os.path.isdir(os.path.join(input_path, p))]
    files = [[os.path.join(input_path, dir, f) for f in os.listdir(os.path.join(input_path, dir)) if f[-3:] == ".gz"][0] for dir in dirs]

    for idx, (name, file) in enumerate(zip(dirs, files)):
        print("\r  - Processing file " + str(idx+1) + "/" + str(len(files)), end="")
        outfile_path = os.path.join(output_path, name + ".fna")

        with gzip.open(file, "rb") as file_in:
            with open(outfile_path, "wb") as file_out:
                shutil.copyfileobj(file_in, file_out)

    print()


def universal_markers(reps, seq_list):
    reps = reps.loc[reps.index.isin(seq_list)]
    mat_reps = reps.to_numpy()

    def universal_markers_inner(mat, threshold = 0.95):
        return np.where(np.count_nonzero(mat == 1, axis=0) / mat.shape[0] >= threshold)[0]

    return universal_markers_inner(mat_reps)


def filter_sequences(pfam_counts, seq_list, markers):
    pfam_counts = pfam_counts.loc[pfam_counts.index.isin(seq_list)]
    pfam_mat = pfam_counts.to_numpy()

    def preestimates_inner(mat, markers_inner):
        submatrix = mat[:, markers_inner]

        completeness = np.sum(np.clip(submatrix, 0, 1), axis=1) / markers_inner.shape[0]
        contamination = np.sum(submatrix - np.clip(submatrix, 0, 1), axis=1) / markers_inner.shape[0]

        return np.array([completeness, contamination], dtype=np.float32).T

    preestimates = preestimates_inner(pfam_mat, markers)
    idx_filtered = np.where(np.logical_and(preestimates[:,0] >= MIN_COMP, preestimates[:,1] <= MAX_CONT))[0]

    return pfam_counts.iloc[idx_filtered]


if __name__ == "__main__":
    main()
