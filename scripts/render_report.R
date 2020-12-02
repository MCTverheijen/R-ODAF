### INSTRUCTIONS

# Path
# Metadata and Contrasts
# Advanced options


# Set custom parameters
params <- list()
params$projectdir <- "~/shared/projects/Template_RNASeq/"
params$species <- "rat"
params$design <- "group" # Must be a single column of the metadata!
params$intgroup <- c("group","dose","sex") # Can be multiple columns of interest (covariates, etc)
params$platform <- "RNA-Seq" # Choose from RNA-Seq or TempO-Seq
# Add parameter here to specify whether pathway analysis should be run
# Add parameter for large sample numbers (or automatically detect >48? 96?)
# Option to limit

# Input file - Rmd
inputFile <- file.path(params$projectdir, "Rmd", "DESeq2_report.rnaseq.Rmd")
# Output file - HTML
outFile <- file.path(params$projectdir, paste0("reports/RNASeq_analysis_", format(Sys.time(), '%d-%m-%Y.%H.%M'),".html"))

rmarkdown::render(input = inputFile,
                  encoding = encoding,
                  output_file = outFile,
                  params = params,
                  envir = new.env())

# For large experiments, loop through contrasts in groups
# - new param - facet_var for "faceting variable"


# Make conda env for ODAF
#conda create -n myenv
#conda install -n myenv scipy=0.15.0