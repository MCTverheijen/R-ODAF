# Set custom parameters
params <- list()
params$projectdir <- "~/shared/projects/Template_RNASeq/"
params$species <- "rat"
params$design <- "group" # Must be a single column of the metadata!
params$intgroup <- c("group","dose","sex") # Can be multiple columns of interest (covariates, etc)

# Input file - Rmd
inputFile <- file.path(params$projectdir, "Rmd", "DESeq2_report.rnaseq.Rmd")
# Output file - HTML
outFile <- file.path(params$projectdir, paste0("reports/RNASeq_analysis_", format(Sys.time(), '%d-%m-%Y.%H.%M'),".html"))

rmarkdown::render(input = inputFile,
                  encoding = encoding,
                  output_file = outFile,
                  params = params,
                  envir = new.env())