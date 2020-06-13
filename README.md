# Sequencing Estimator
## Rhett M. Rautsaw

This Shiny application is designed to help you choose the right sequencing platform and estimate costs.

It is designed to estimate sequencing needs and costs for:

- Whole Genome Sequencing (WGS)
- Transcriptomics (RNA-seq)
- Whole Genome Bisulfite Sequencing (WGBS)
- Assay for Transposase Accessible Chromatin (ATAC-Seq)
- Double-digest Restriction-site Associated DNA Sequencing (ddRAD-Seq)
- **However**, given a desired coverage or number of reads (in millions), the different areas could be used for **any** type of sequencing.

Most of the information for sequencing platorms was obtained from [Illumina](https://www.illumina.com/systems/sequencing-platforms.html). Illumina, PacBio, and Nanopore are constantly updating their platforms and I may not be able to keep up. Therefore, the option is provided to download an example database and update it yourself to re-upload.

The recommended number of reads and coverage for each type of sequencing are default input for all options (*e.g.* 20x for WGBS, 5 M reads for ddRAD-Seq). 

## Running the Application


