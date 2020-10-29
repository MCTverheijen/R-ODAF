######################################################
### DEseq2 for normalization and DE identification ###
######################################################

library(DESeq2)
#require("pheatmap")
#require("ggplot2")
require("DESeq2")
require("edgeR")
#library("lattice")


###################################################################################
###################################################################################
# PARAMETERS TO SET MANUALLY                            		 

# Set file locations				
	sampledir <- "/ngs-data-2/data/CEFIC/CEFIC_scripts/DESeq2_Github/TestData/"
	outputdir <- paste(sampledir, "Output/", sep="")
	if(!dir.exists(outputdir)) {dir.create(outputdir)}
# Names of files to load
	SampleDataFile <- "sampleData.csv" #This comma delimited file contains the merged RSEM.genes.results files
	SampleKeyFile <- "samplekeyTEST.csv" #This comma delimited file contains at least 2 columns: NAME (sample names identical to the column names of sampleData) and Compound (needs to identify to which group the sample belongs -> ExperimentalGroup & ControlGroup)

# Specify which groups need to be compared 
	Samp4compare<- c("PIR") # Experimental group, needs to correspond with the SampleKeyFile
	Cont4compare<- c("Vehicle") # Control group, needs to correspond with the SampleKeyFile
	DESIGN<- "Compound"	#Column name samplekeyTEST.csv which defines the groups to be compared

# Set analysis ID. This ID will be used as prefix for the output files
	analysisID <-"R-ODAF_test_PIR"
# Specify used platform/technology for data generation:
	Platform <- "RNA-Seq" # Specify "RNA-Seq" or "TempO-seq"





###################################################################################
#DEFINE FUNCTIONS
###################################################################################
plot.barplots<-function(samples,b) {
	color <- NULL
	for (h in 1:ncol(norm_data)){
		if (substring(colnames(norm_data)[h], 1, 3) == substring(Samp4compare, 1, 3)) { color <- c(color, "red3") } else { color <- c(color, "darkgrey")}
	}
	fileNamePlot <- paste0(b, row.names(samples)[b], ".png")
	pseudoTitle <- paste0(row.names(samples)[b], "_pAdj:", samples[b,"padj"])
	      
	png(file=paste(fileNamePlot, sep="/"), width=1200, height=700, pointsize=20)
		par(mar=c(8,4,3,1))
	    	barplot(as.numeric(norm_data[row.names(samples)[b],]), las=2, col=color, main=pseudoTitle, cex.names=0.5,  cex.axis=0.8, names.arg=colnames(norm_data)) 
	dev.off()
} #plot.barplots function done

###################################################################################
draw.barplots<-function(samples, top_bottom, NUM){
	if (nrow(samples) == 0) {
		#print("no genes to plot") 
	} else { 
	if (top_bottom == "top") {
		#print(paste0("drawing Top ", NUM, " plots"))
		if (nrow(samples) <= NUM) { 
			for (b in 1:nrow(samples)) {plot.barplots(samples,b)}
		}

		if (nrow(samples) > NUM) { 
			for (b in 1:NUM) {plot.barplots(samples,b)}	
		}	
	}

	if (top_bottom == "bottom") {
		#print(paste0("drawing Bottom", NUM, " plots"))
		if (nrow(samples) <= NUM) { 
			for (b in 1:nrow(samples)) {plot.barplots(samples,b)}
		}
		if (nrow(samples) > NUM) { 
			for (b in ((nrow(samples)-NUM+1):nrow(samples))) {plot.barplots(DEsamples,b)}
		}
	}}
} #draw.barplots function done

###################################################################################
###################################################################################

#Set parameters according to platform
if (Platform=="RNA-Seq"){
	MinCount<- 1
	pAdjValue<- 0.01 
} else if (Platform=="TempO-seq") {
	MinCount<- 0.5
	pAdjValue<- 0.05 
} else { print("Platform/technology not recognized") }



# Load input files 
setwd(sampledir)
sampleData <- read.delim(SampleDataFile, sep=",", stringsAsFactors=FALSE, header=TRUE,  quote="\"", row.names=1)
DESeqDesign <- read.delim(SampleKeyFile, stringsAsFactors=FALSE, sep=",", header=TRUE,  quote="\"", row.names="NAME")

NORM_TYPE<-paste0(analysisID, "_DESeq2_", Platform)
print(NORM_TYPE)

plotdir<- paste(outputdir, "/plots/", sep="")
if(!dir.exists(plotdir)) {dir.create(plotdir)}
barplot.dir<- paste(plotdir, "/barplot_genes/", sep="")
if(!dir.exists(barplot.dir)) {dir.create(barplot.dir)}

# First data clean-up: replace NA & remove samples with total readcount < threshold
sampleData[ is.na(sampleData) ] <- 0 
sampleData<- sampleData[,(colSums(sampleData)>1000000)]

##########
# DESeq2 #
##########

for (x in 1:length(Samp4compare)){	## for all comparisons to be done	
	condition1<- Cont4compare[x]	    		
	condition2<- Samp4compare[x]  

	DE_Design <- matrix(data=NA, ncol=2)
	DE_Design <- DESeqDesign [c(grep(condition1,DESeqDesign[,DESIGN[x]]), grep(condition2,DESeqDesign[,DESIGN[x]])),]
	samples <- sampleData[, rownames(DE_Design) ]

	###########
	print(paste(condition2, " vs ", condition1, ":", NORM_TYPE))		

	colnames(samples)<-NULL
	if (DESIGN[x] == "Compound") {
		dds <- DESeqDataSetFromMatrix(countData = round(samples), colData = as.data.frame(DE_Design), design = ~ Compound)
	} else {print("Setting of DESIGN failed. Please check code: dds <- DESeqDataSetFromMatrix(countData = round(samples), colData = as.data.frame(DE_Design), design = ~ Compound)}.")}

	print("Wait... (dds step executing)")				
	dds <- DESeq(dds, quiet=TRUE)
	#Filter low readcounts 
	print(paste0("Filtering genes: 75% of at least 1 group need to be above ", MinCount, " CPM"))
	print("AND")
	print("Detecting spurious spikes: Max-Median > Sum/(Rep+1)" )
		SampPerGroup<-table(DE_Design[,DESIGN])
		idx<-FlagSpike<-NameRows<-NULL
		Counts<-counts(dds, normalized=TRUE)
		CPMdds<-cpm(counts(dds, normalized=TRUE))
		for (gene in 1:nrow(dds)) {
			GroupsPass<-checkSpike<-NULL
			for (group in 1:length(SampPerGroup)) { #test if group passes
				sampleCols<-grep(dimnames(SampPerGroup)[[1]][group],DE_Design[,DESIGN])
				Check<-sum(CPMdds[gene,sampleCols] >= MinCount)>= 0.75*SampPerGroup[group]
				GroupsPass<-c(GroupsPass, Check)
				if (Check == FALSE) {checkSpike<- c(checkSpike, Check)} else {
					checkSpike<-c(checkSpike, ((max(Counts[gene,sampleCols])-median(Counts[gene,sampleCols])) >= (sum(Counts[gene,sampleCols])/(SampPerGroup[group]+1))))
				}
			}
			idx <- c(idx, as.logical(sum(GroupsPass)))
			if (sum(checkSpike) >=1) {
				FlagSpike<-rbind(FlagSpike, Counts[gene,])
				NameRows<<-c(NameRows, row.names(Counts)[gene])
				row.names(FlagSpike)<-NameRows 
			}
		}		

	print("Obtaining the DESeq2 results")
	res <- results(dds[idx], contrast=c(DESIGN[x], condition2, condition1), pAdjustMethod= 'fdr')   
	setwd(outputdir)
	FileName<-paste(NORM_TYPE, condition2,"vs",condition1, "FDR", pAdjValue, sep="_")		
	#Save output tables		
	norm_data <<- counts(dds[idx],normalized=TRUE) 
	write.table(norm_data,file=paste0(FileName, "_Norm_Data.txt"), sep="\t", quote=FALSE)
	write.table(FlagSpike,file=paste0(FileName, "_FlaggedSpikes.txt"), sep="\t", quote=FALSE)
	DEsamples <<- subset(res,res$padj < pAdjValue)	
	write.table(DEsamples,file=paste0(FileName,"_DEG_table.txt"), sep="\t", quote=FALSE)
	DEspikes<<- DEsamples[rownames(DEsamples)%in%NameRows,]	
	write.table(DEspikes,file=paste0(FileName,"_DEspikes_table.txt"), sep="\t",quote=FALSE)

	print("creating Read count Plots")
	# top DEGs
	plotdir<- paste(outputdir, "/plots/", sep="")
	if(!dir.exists(plotdir)) {dir.create(plotdir)}
	barplot.dir<- paste(plotdir, "/barplot_genes/", sep="")
	if(!dir.exists(barplot.dir)) {dir.create(barplot.dir)}

	TOPbarplot.dir<- paste(barplot.dir, "Top_DEGs/", sep="")
	if(!dir.exists(TOPbarplot.dir)) {dir.create(TOPbarplot.dir)}
	setwd(TOPbarplot.dir)
	draw.barplots(DEsamples, "top", 20) #(DEsamples, top_bottom, NUM)
	print("Top 20 DEG plots done")

	# Spurious spikes
	SPIKEbarplot.dir<- paste(barplot.dir, "DE_Spurious_spikes/", sep="")
	if(!dir.exists(SPIKEbarplot.dir)) {dir.create(SPIKEbarplot.dir)}
	setwd(SPIKEbarplot.dir)
	draw.barplots(DEspikes, "top", nrow(DEspikes)) #(DEsamples, top_bottom, NUM)
	print("All DE_Spurious_spike plots done")


	print("DESeq2 Done")
}
print("END of script. Have a nice day!")




