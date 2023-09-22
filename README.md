# CoCoPyE Database

This is the database repository for CoCoPyE. It is mainly used as a download source. If you are looking for the tool itself, see https://github.com/gobics/cocopye.

## Database Structure

CoCoPyE can be used with Pfam versions 24 and 28. It has to use a database that was build with the same version. They can be found in the corresponding folders (`24` and `28`). `version.txt`contains the current version (has to be incremented manually when a new release is published). `info.txt` contains information about the dataset source.

## Database Build Process

The following section explains how to build the database for some Pfam version (which results in the content of folder `24`/`28`).

### Input

You need a folder containing the sequences you want to use for the database (depending on how you choose them, these can be the same or different for both Pfam versions), one FASTA file for each sequence. FASTA headers in the files are ignored. You also need a csv-file with metadata, one row for each input sequence. The required columns are `sequence`, `superkingdom`, `phylum`, `class`, `order`, `family`, `genus` and `species` where `sequence` is the name of the corresponding FASTA file and the other columns are taxonomy data. All cells except `sequence` are allowed to be empty. However, the quality of the taxonomy prediction is directly related to the quality of the metadata. Therefore it is recommended to avoid empty cells if possible.

### Building the database

The build process of the database is part of the CoCoPyE CLI. Just run `cocopye database -i <folder with FASTA files> -m <metadata-file> -o <output folder>`. If you want to build with Pfam 24, use `cocopye --pfam24 database -i <folder with FASTA files> -m <metadata-file> -o <output folder>`.

### Additional files

The build command automatically creates most required files. However, you still have to manually add the files for the machine learning method (most likely a neural network). The files have to be named `model_comp.pickle` for the completeness predictor and `model_cont.pickle` for contamination. If you want to make a new release, you should also increment the version number in `version.txt` to match the release version.