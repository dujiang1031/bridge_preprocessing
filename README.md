# HSCmanuscript
# python scripts for analysis in manuscript "Variations in Gene Expression Distinguish in vivo Cellular Activity Levels "


1. construct_molecular-bridges.py

  the script to step-wise construct the bridge connecting cell-barcode (CBC) and viral clonal tracking barcode (TBC)


2. 00-external

  folder containing supporting materials for data processing, including the TBC library IDs, viral backbone sequences, multiplex information, etc


3. required additional input:

  matrix data from clonal tracking;
  fastq files from molecular bridge sequencing by PacBio;
  H5 file from single cell RNAseq processed by cellranger.
