####################################################################################
### Sequencing R-ODAF, Omics Data Analysis Framework for Regulatory application  ###
####################################################################################

####################################################
#### Settings which need to be adapted by user #####
####################################################

# ! NOTE ! # The docker is mounted on a specific folder of your system. This part of the file path is replaced by /data/.
# Example: Mounted on        /home/username/project/
#          Samples are in    /home/username/project/Samples/
#          Script should say /data/Samples/

#specify the directory for the output
OUTPUT_DIR="/data/test/output/"
#specify location of input fastq files. ALL FILES IN THE FOLDER WILL BE PROCESSED 
RAW_SAMPLE_DIR="/data/test/"
# specify extention of input files (".fastq" or ".fastq.gz") 
SUFFIX_INPUTFILES='.fastq.gz' 
#specify the sequencing mode used to obtain the data
SEQMODE="single" #specify "paired" or "single" end mode
# specify the read suffix (e.g. "_R1_001")
PAIRED_END_SUFFIX_FORWARD="_R1_001"
# *IF* paired end mode was used, specify the reverse suffix as well (e.g. "_R2")
PAIRED_END_SUFFIX_REVERSE="_R2"
#
# Choose the main organism for genome alignment (e.g "Rat_6.0.97"). {NOTE: This ID is a lable specific for this script and is made for the user to identify which eneome wersion was used. It can contain any text}.
ORGANISM_GENOME_ID="Rat_6.0.97" 
# PATH/Directory in which the genomefiles are located
GENOME_FILES_DIR="/data/GENOME/"
# Filename of genome fasta file (without path)
GENOME_FILE_NAME="Rattus_norvegicus.Rnor_6.0.dna.toplevel.fa"
# Filename of GTF file (without path)
GTF_FILE_NAME="Rattus_norvegicus.Rnor_6.0.97.gtf"
# Whether the genome indexing has already been done. When "Yes" is specified, the indexing will be skipped. If "No" The index will be made
GENOME_INDEX_ALREADY_AVAILABLE="Yes" #Specify "Yes" or "No"
#Specify whether you are working with a large (=human) genome. Specify "Yes" when working with human or "No"
LARGE_GENOME="No"
#
# System parameters
#Specify amount of CPUs to use for alignment step (advised=20 or 30)
CPU_FOR_ALIGNMENT=20 
#Specify amount of CPUs to use (advised: 6 or higher)
CPU_FOR_OTHER=6 

### No other input required ###

#################################
# Running the sequencing R-ODAF #
#################################

declare BASEDIR=${OUTPUT_DIR}   #specify workdirectory
#declare SEQMODE="single"   #specify "paired" or "single" end mode
#declare RAW_SAMPLE_DIR=${BASEDIR}   #specify location of fastq files
declare SUFFIX_IN=${SUFFIX_INPUTFILES} # specify extention of input files (".fastq" or ".fastq.gz")

#For paired end mode: specify the SUFFIXES used for read1 and read2 
declare PAIR1=${PAIRED_END_SUFFIX_FORWARD}
declare PAIR2=${PAIRED_END_SUFFIX_REVERSE}

#Specify reference genome files
declare GENOMEDIR=${GENOME_FILES_DIR}       #location of genomefiles
declare GENOME="${GENOMEDIR}${GENOME_FILE_NAME}"  #genome fasta file
declare	GTF="${GENOMEDIR}${GTF_FILE_NAME}"   #genome GTF file
declare GenomeID=${ORGANISM_GENOME_ID}  #Specify the genome name (e.g. Species+GenomeVersion or Species+DownloadDate) to prevent overwriting other indexed genomes

declare LargeGenome="No"     #Specify "Yes" or "No" to indicate if Human=Large Genome is used
GenomeIndexDone=${GENOME_INDEX_ALREADY_AVAILABLE}     #Specify "Yes" or "No" to indicate if GenomeIndexing has already been done (if no is specified, the index will be made)

declare CPUs=${CPU_FOR_OTHER}   #Specify amount of CPUs to use
declare CPUs_align=${CPU_FOR_ALIGNMENT}   #Specify amount of CPUs to use for alignment step

#######################################################################
### Defining parameters for script to run (no user input nessesary) ### 
#######################################################################

declare SOURCEDIR=${RAW_SAMPLE_DIR}
declare TRIMM_DIR="${BASEDIR}/Trimmed_reads/"
declare OUTPUTDIR=${TRIMM_DIR}
declare QC_DIR_fastp="${OUTPUTDIR}/fastpQCoutput/"
declare QC_DIR_multiQC="${OUTPUTDIR}/MultiQC/"
declare align_DIR="${OUTPUTDIR}STAR"
declare Quant_DIR="${OUTPUTDIR}RSEM/"
declare RSEM_GENOMEDIR="${GENOMEDIR}RSEM/"

declare SUFFIX1=${SUFFIX_IN} 
declare SUFFIX_out="_trimmed${SUFFIX_IN}"

mkdir -p ${TRIMM_DIR}
mkdir -p ${QC_DIR_fastp}
mkdir -p ${QC_DIR_multiQC}
mkdir -p ${align_DIR}
mkdir -p ${Quant_DIR}
mkdir -p ${RSEM_GENOMEDIR}

###################################################################################################
### Trimming raw reads : Fastp ###
##################################

source activate fastp
# trimming single end reads  
if [ ${SEQMODE} == "single" ]; then
	declare FILES1="${SOURCEDIR}*${SUFFIX1}"
for FILENAME in ${FILES1[@]}; do
	echo -e "[TRIMMING] fastp: [${FILENAME:${#SOURCEDIR}:-${#SUFFIX1}}]" 
	#Single end
	fastp --in1 ${FILENAME} --out1 ${TRIMM_DIR}${FILENAME:${#SOURCEDIR}:-${#SUFFIX1}}${SUFFIX_out} --json "${QC_DIR_fastp}${FILENAME:${#SOURCEDIR}:-${#SUFFIX1}}_fastp.json" --html ${QC_DIR_fastp}"${FILENAME:${#SOURCEDIR}:-${#SUFFIX1}}fastp.html" --cut_front --cut_front_window_size 1 --cut_front_mean_quality 3 --cut_tail --cut_tail_window_size 1 --cut_tail_mean_quality 3 --cut_right --cut_right_window_size 4 --cut_right_mean_quality 15 --length_required 36
done; fi

# trimming paired end reads  
if [ ${SEQMODE} == "paired" ]; then
	declare FILES1="${SOURCEDIR}*${PAIR1}${SUFFIX1}";
for FILENAME in ${FILES1[@]}; do
	READ1=${FILENAME}	
	READ2=${FILENAME:0:-${#SUFFIX1}-${#PAIR1}}${PAIR2}${SUFFIX1}
	echo -e "[TRIMMING] fastp: [${READ1:${#SOURCEDIR}:-${#PAIR1}-${#SUFFIX1}}]" 
	fastp --in1 ${READ1} --in2 ${READ2} --out1 ${TRIMM_DIR}${READ1:${#SOURCEDIR}:-${#SUFFIX1}}${SUFFIX_out} --out2 ${TRIMM_DIR}${READ2:${#SOURCEDIR}:-${#SUFFIX1}}${SUFFIX_out} --json "${QC_DIR_fastp}${READ1:${#SOURCEDIR}:-${#SUFFIX1}-${#PAIR1}}PE_fastp.json" --html "${QC_DIR_fastp}${READ1:${#SOURCEDIR}:-${#SUFFIX1}-${#PAIR1}}PE_fastp.html" --cut_front --cut_front_window_size 1 --cut_front_mean_quality 3 --cut_tail --cut_tail_window_size 1 --cut_tail_mean_quality 3 --cut_right --cut_right_window_size 4 --cut_right_mean_quality 15 --length_required 36 
done; fi
conda deactivate

################################
#INFORMATION ON TRIMMING PROCESS
#--cut_front --cut_front_window_size 1 --cut_front_mean_quality 3 == Trimmomatic  "LEADING:3"
#--cut_tail --cut_tail_window_size 1 --cut_tail_mean_quality 3 == Trimmomatic  "TRAILING:3"
#--cut_right --cut_right_window_size 4 --cut_right_mean_quality 15 == Trimmomatic  "SLIDINGWINDOW:4:15"

#### If additional trimming is needed (see multiQCreport): 
# add to R1: --trim_front1 {amount_bases} and/or --trim_tail1 {amount_bases}
# add to R2: --trim_front2 {amount_bases} and/or --trim_tail2 {amount_bases}


###################################################################################################
### Quality control raw reads: Fastp + MultiQC report ###
################################################################

source activate multiqc
# Running multiQC on fastp-output
multiqc ${QC_DIR_fastp} --filename MultiQC_Report.html --outdir ${QC_DIR_multiQC}
conda deactivate

####################################################
### Alignment of reads (paired end & single end) ###
####################################################

source activate star
#Indexing reference genome for STAR
cd ${GENOMEDIR}
if [ ${GenomeIndexDone} == "No" ]; then
	if [ ${LargeGenome} == "Yes" ]; then
	echo "Indexing Large(=Human) genome"
	# additional parameters to save RAM usage: --genomeSAsparseD 1 --genomeChrBinNbits 15
	STAR --runMode genomeGenerate --genomeDir ${GENOMEDIR} --genomeFastaFiles ${GENOME} --sjdbGTFfile ${GTF} --sjdbOverhang 99 --runThreadN ${CPUs} --genomeSAsparseD 1 --genomeChrBinNbits 15
	else
	echo "Indexing Small(=non-human) genome"
	STAR --runMode genomeGenerate --genomeDir ${GENOMEDIR} --genomeFastaFiles ${GENOME} --sjdbGTFfile ${GTF} --sjdbOverhang 99 --runThreadN ${CPUs}
	fi
fi
cd ${SOURCEDIR}

# Aligning reads single end
if [ ${SEQMODE} == "single" ]; then
	declare FILES1=${OUTPUTDIR}*${PAIR1}*${SUFFIX1}		
for FILENAME in ${FILES1[@]}; do
	READ1=${FILENAME:0:-${#SUFFIX1}}	
	echo -e "[ALIGNING] STAR : [${READ1:${#OUTPUTDIR}:-${#PAIR1}}]" 
	if [ ${SUFFIX1} == ".fastq.gz" ]; then
	gunzip ${READ1}${SUFFIX1}
	fi
 	STAR --runThreadN ${CPUs_align} --genomeDir ${GENOMEDIR} --readFilesIn "${READ1}.fastq" --quantMode TranscriptomeSAM --outFileNamePrefix ${align_DIR}${READ1:${#OUTPUTDIR}:-(${#PAIR1}+8)}
done;fi

# Aligning reads paired end
if [ ${SEQMODE} == "paired" ]; then
	declare FILES1=${OUTPUTDIR}*${PAIR1}*${SUFFIX1}		
for FILENAME in ${FILES1[@]}; do
	READ1=${FILENAME:0:-${#SUFFIX1}}	
	READ2=${FILENAME:0:-${#SUFFIX1}-(${#PAIR1}+8)}${PAIR2}"_trimmed"
	echo -e "[ALIGNING] STAR : [${READ1:${#OUTPUTDIR}:-(${#PAIR1}+8)}]" 
	if [ ${SUFFIX1} == ".fastq.gz" ]; then
		gunzip ${READ1}${SUFFIX1}
		gunzip ${READ2}${SUFFIX1}
	fi
 	STAR --runThreadN ${CPUs_align} --genomeDir ${GENOMEDIR} --readFilesIn "${READ1}.fastq" "${READ2}.fastq" --quantMode TranscriptomeSAM --outFileNamePrefix ${align_DIR}${READ1:${#OUTPUTDIR}:-(${#PAIR1}+8)}
done;fi
conda deactivate

#######################
# QUANTIFICATION RSEM #
#######################

source activate rsem
#Indexing reference genome for RSEM
cd ${GENOMEDIR}
if [ ${GenomeIndexDone} == "No" ]; then 
	rsem-prepare-reference --gtf ${GTF} ${GENOME} ${RSEM_GENOMEDIR}${GenomeID}
	GenomeIndexDone="Yes"
fi

cd ${SOURCEDIR}
# Quantify reads single end
if [ ${SEQMODE} == "single" ]; then
	declare FILES1=${OUTPUTDIR}*${PAIR1}*${SUFFIX1}		
for FILENAME in ${FILES1[@]}; do
	READ1=${FILENAME:0:-${#SUFFIX1}}	
	echo -e "[QUANTIFYING] RSEM : [${READ1:${#OUTPUTDIR}:-(${#PAIR1}+8)}]" 
	rsem-calculate-expression -p ${CPUs} --bam "${align_DIR}${READ1:${#OUTPUTDIR}:-(${#PAIR1}+8)}Aligned.toTranscriptome.out.bam" --no-bam-output ${RSEM_GENOMEDIR}${GenomeID} ${Quant_DIR}${READ1:${#OUTPUTDIR}:-(${#PAIR1}+8)}
done;fi

# Quantify reads paired end
if [ ${SEQMODE} == "paired" ]; then
	declare FILES1=${OUTPUTDIR}*${PAIR1}*${SUFFIX1}		
for FILENAME in ${FILES1[@]}; do
	READ1=${FILENAME:0:-${#SUFFIX1}}	
	READ2=${FILENAME:0:-${#SUFFIX1}-(${#PAIR1}+8)}${PAIR2}"_trimmed"
	echo -e "[QUANTIFYING] RSEM : [${READ1:${#OUTPUTDIR}:-(${#PAIR1}+8)}]" 
	rsem-calculate-expression -p ${CPUs} --paired-end --bam "${align_DIR}${READ1:${#OUTPUTDIR}:-(${#PAIR1}+8)}Aligned.toTranscriptome.out.bam" --no-bam-output ${RSEM_GENOMEDIR}${GenomeID} ${Quant_DIR}${READ1:${#OUTPUTDIR}:-(${#PAIR1}+8)}
done;fi
conda deactivate

#Rezip unzipped files
if [ "${SUFFIX1}" == ".fastq.gz" ]; then
	rm ${OUTPUTDIR}*"trimmed.fastq"; fi
echo "Pre-processing of data complete"

