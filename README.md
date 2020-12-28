# HSCmanuscript
# python scripts for analysis in manuscript "Quantifying in vivo associations between gene expression and stem cell function at the single cell level"


1. construct_molecular-bridges.py

  the script to step-wise construct the bridge connecting cell-barcode (CBC) and viral clonal tracking barcode (TBC).
  To use the code, follow these steps:
    1) download the ".py" code file and the "00-external" folder into the same directory
    2) prepare the additional input data (see below)
    3) open the ".py" code in a Python interpreter (Spyder recommended)
    4) uncomment the session and run


2. 00-external

  folder containing supporting materials for data processing, including the TBC library IDs, viral backbone sequences, multiplex information, etc


3. required additional input:

  matrix data from clonal tracking;
  fastq files from molecular bridge sequencing by PacBio;
  H5 file from single cell RNAseq processed by cellranger.
