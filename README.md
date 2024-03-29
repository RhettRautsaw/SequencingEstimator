# Sequencing Estimator
<img align="right" src="www/SequencingEstimator.png" width=150>

## Rhett M. Rautsaw

[https://RhettRautsaw.app/shiny/SequencingEstimator/](https://RhettRautsaw.app/shiny/SequencingEstimator/)

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

There are many ways to download and run it. 

It is hosted freely on my Digital Ocean Droplet: 
[https://reptilerhett.app/shiny/SequencingEstimator/](https://reptilerhett.app/shiny/SequencingEstimator/)

This app can also be run through R:

```R
library(shiny)

# Easiest way is to use runGitHub
runGitHub("SequencingEstimator", "RhettRautsaw")

# Run a tar or zip file directly
runUrl("https://github.com/RhettRautsaw/SequencingEstimator/archive/master.tar.gz")
runUrl("https://github.com/RhettRautsaw/SequencingEstimator/archive/master.zip")
```

Or you can clone the git repository, then use `runApp()`:

```R
# First clone the repository with git. If you have cloned it into
# ~/SequencingEstimator, first go to that directory, then use runApp().
setwd("~/SequencingEstimator")
runApp()
```


To run a Shiny app from a subdirectory in the repo or zip file, you can use the `subdir` argument. This repository happens to contain another copy of the app in `inst/shinyapp/`.

```R
runGitHub("SequencingEstimator", "RhettRautsaw", subdir = "inst/shinyapp/")

runUrl("https://github.com/RhettRautsaw/SequencingEstimator/archive/master.tar.gz",
  subdir = "inst/shinyapp/")
```
