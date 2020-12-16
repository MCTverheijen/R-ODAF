#!/usr/bin/R
# Custom parameters for the report
library(tidyverse)
params <- list()
params$projectdir <- here::here() # If you open this as an R project, it should find your root directory
params$project_name <- "Bisphenols"
params$species <- "human"     # one of human, mouse, rat, hamster
params$design <- "group"    # single experimental group of interest
params$intgroup <- c("group","dose","pool","batch") # Can be multiple columns of interest (covariates, etc)
params$flag <- TRUE         # runs all analysis by default when set to TRUE
params$platform <- "TempO-Seq" # Choose from RNA-Seq or TempO-Seq
params$group_facet <- "chemical" # NULL # If you have many different experimental groups, you may subset the report by specifying a column in the metadata to filter groups, and then setting the group of interest in group_filter
params$exclude_samples <- NA  # Optionally, a vector of sample names to exclude from the analysis
params$exclude_groups <-  c("cells, no treatment",
                            "No Cells; Brain +ve control",
                            "No cells; media only",
                            "No cells; Negative control",
                            "No cells; Universal human +ve control") # Optionally, a vector of groups to exclude from the analysis. By default this is assumed to be in the column specified by params$design.
params$include_only_column <- NA # Restrict analysis to group(s) in the column listed here based on params$include_only_group.
params$include_only_group <- NA # Restrict analysis to this/these group(s) within the column listed in params$include_only_column
params$use_cached_RData <- FALSE # If possible, load the saved RData for dds object and gene IDs
params$cpus <- 39 # Set to a lower number (e.g., 2 to 4) if you aren't working in a server environment
params$run_pathway_analysis <- TRUE

skip_extra <- "DMSO 0.1%"

# Input file - Rmd
inputFile <- file.path(params$projectdir, "Rmd", "DESeq2_report.rnaseq.Rmd")

if(is.na(params$group_facet)){
  # Output file - HTML
  outFile <- file.path(params$projectdir, paste0("reports/RNASeq_analysis_", params$project_name, "_",
                                                 format(Sys.time(),'%d-%m-%Y.%H.%M'),".html"))
  rmarkdown::render(input = inputFile,
                    encoding = encoding,
                    output_file = outFile,
                    params = params,
                    envir = new.env())
} else {
  SampleKeyFile <- file.path(params$projectdir,"metadata/metadata.txt")
  DESeqDesign <- read.delim(SampleKeyFile,
                            stringsAsFactors=FALSE,
                            sep="\t",
                            header=TRUE,
                            quote="\"",
                            row.names=1) # Pick column that is used in ID; might be more appropriate to change this!
  DESeqDesign$original_names <- rownames(DESeqDesign)
  facets <- DESeqDesign %>% filter(!(!!sym(params$group_facet)) %in% c(params$exclude_groups, skip_extra)) %>% 
    pull(params$group_facet) %>% 
    unique() # remove params$exclude_groups
  message(paste0("Making multiple reports based on ", params$group_facet ,"..."))
  for(i in facets[12:length(facets)]){
    message(paste0("Building report for ", i, "..."))
    params$group_filter <- i
    outFile <- file.path(params$projectdir, paste0("reports/RNASeq_analysis_", params$project_name, "_",
                                                   i,
                                                   "_",
                                                   format(Sys.time(),'%d-%m-%Y.%H.%M'),".html"))
    rmarkdown::render(input = inputFile,
                      encoding = encoding,
                      output_file = outFile,
                      params = params,
                      envir = new.env())
  }
}
