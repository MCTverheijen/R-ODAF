# Custom parameters for the report
params <- list()
params$projectdir <- "~/shared/projects/Template_RNASeq/"
params$species <- "rat"     # one of human, mouse, rat, hamster
params$design <- "group"    # single experimental group of interest
params$intgroup <- c("group","dose","sex") # Can be multiple columns of interest (covariates, etc)
params$flag <- TRUE         # runs all analysis by default when set to TRUE
params$platform <- "RNA-Seq" # Choose from RNA-Seq or TempO-Seq
params$group_facet <- "sex"  # If you have many different experimental groups, you may subset the report by specifying a column in the metadata to filter groups, and then setting the group of interest in group_filter
params$group_filter <- "m"           # Which group will this report be done on?
params$exclude_samples <- c("SRR5890440","SRR5890376")  # Optionally, a vector of sample names to exclude from the analysis
params$exclude_groups <- c("f2","m2") # Optionally, a vector of groups to exclude from the analysis. By default this is assumed to be in the column specified by params$design.

# Input file - Rmd
inputFile <- file.path(params$projectdir, "Rmd", "DESeq2_report.rnaseq.Rmd")

# Output file - HTML
outFile <- file.path(params$projectdir, paste0("reports/RNASeq_analysis_", format(Sys.time(), '%d-%m-%Y.%H.%M'),".html"))

if(is.null(params$group_facet)){
  rmarkdown::render(input = inputFile,
                    encoding = encoding,
                    output_file = outFile,
                    params = params,
                    envir = new.env())
} else {
  message(paste0("Making multiple reports based on ", params$group_facet,"..."))
  for
}
