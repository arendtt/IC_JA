# ENIGMA - CNVs

## CNV detection and QS

### **1. Organize your data and identify appropriate files for CNV calling**
Create a folder called “ENIGMA-CNV_Analysis/” where you wish your output and software for the ENIGMA-CNV protocol to be. Henceforward referred to as “the Analysis-folder”.

OBS – note that if you have data from more than one type of chip, please run the protocol separately for each different chip set
Please place all generated files in the Analysis-folder.

a. Genetic information files
CNVs are called based on data from genotyping chips:

### **a.1 Illumina - IlluminaFinalReport.txt/LRR-BAF-files**

-the **Illumina Final report file** of your Dataset/samples with, as a minimum:
 
[Header]
BSGT Version    3.2.23
Processing Date 10/31/2008 11:42 AM
Content         sample.bpm
Num SNPs        45707
Total SNPs      45707
Num Samples     48
Total Samples   192
[Data]
SNP Name       Sample ID       B Allele Freq   Log R Ratio
rs1000000  KS2231000715    1.0000  -0.0558
rs1000002  KS2231000715    1.0000  -0.0422
Etc…
OBS – file may not contain the first 8 lines; a previously extracted Illumina Final Report with (as minimum) the mentioned columns is also perfect.

Please place the IlluminaFinalReport.txt-file in the Analysis-folder.

- Alternatively (to the IlluminaFinalReport.txt-file), **LRR-BAF files**. Please note that if these files are present in several subfolders, the script might run into issues – please contact the helpdesk for help if this is the case.

Format of (tab-delimited) LRR-BAF file [Name: SubjID1]:
Name    SubjID1.Log R Ratio           SubjID1.B Allele Freq
rs1000000   -0.0038  0.0161
rs1000002   0.0073   0.9943
rs10000023  -0.0307  0.0026
etc…

#### * Iremos utilizar este LRR-BAF files, mas estará em subfolders 

### **SNP Position File** 
For generation of the PFB-file (population frequency file, see later), you need to provide a SNP-position file.
The SNP-position-file is a tab-delimited file with the positions of the SNPs on your chip, containing at least these columns:
Name    Chromosome       Position
rs1000000   12  126890980
rs1000002   3   183635768
rs10000023  4   95733906

### **b. Cohort-generated files**
### **Sex - SexFile.txt**
-please produce a tab-delimited file specifying SampleID (as in the lluminaFinalReport-file) and sex in the format:
SampleID    sex
SampleID1   female
SampleID2   male
[please note that sex MUST be in the format "female" or "male"; if missing sex-info, just leave the individual out of the list]. This sex-info is also needed for the covariate- file below.
If sex is missing for some individuals, please note the number down in the CNV calling script.
Place the file in the Analysis-folder.

#### * Fazer no RStudio

### **Individuals to remove - RemoveFile.txt**
Please produce a text-file with the individuals you want to remove  (if any) in the analysis.
Each row has the SampleID [and ONLY this ID] specified as in the IlluminaReportFile
SampleID
SampleID1
SampleID2
The RemoveFile.txt should contain individuals that you know are sample mixups or subjects that for instance withdrew from the study so that these are kept out of the analysis. Please note that these are the ONLY initial exclusion measures we wish you to apply as some CNVs are so low in frequency that we want to be as inclusive as possible.7
We request that you keep both duplicates and related individuals (and indicate these in the duplication/relatedness file below) since these will be helpful in QC as well as in analysis. 
Place the file in the Analysis-folder

#### * Não se aplica 

### **Individuals to keep - KeepFile.txt**
Needed if you have files without imaging in your IlluminaFinalReport-file (or for others reasons just want to keep some individuals in the analysis)
-produce a text-file with the individuals you want to keep in the analysis Each row has the SampleID [and ONLY this ID] specified as in the IlluminaReportFile:
SampleID
SampleID1
SampleID2
Place the file in the Analysis-folder

### **Duplication/relatedness DupsRelatives.txt**
Please produce a tab-delimited list of duplicate individuals and related individuals:
SampleID1   SampleID2   PI_HAT  Relation
SampleID10  SampleID24  0.5 Mother_child

-where SampleID1 and SampleID2 are the related individuals, PI_HAT the proportion of the genomic variation shared IBD and Relation (if known or suspected - e.g. duplicate, sibling, parent/child, uncle/nephew, grandparent/grandchild, halfsibling, cousin if known). If Pihat is unknown, keep Pihat as “_“, but instead indicate Relation. Please only include PI_HAT>0.125. 
Pihat can be calculated with the plink –genome command 
(https://zzz.bwh.harvard.edu/plink/ibdibs.shtml). 
Place the file in the Analysis-folder.

#### * Como só os probandos tem neuroimagem, não iremos utilizar os pais, portanto não iremos fazer esse arquivo

### **2. Download necessary software, files and scripts**
Please place all scripts and files downloaded in the Analysis-folder.

### **a. Container software (independent of dataset)**
To minimize the impact of differences in computer hardware and software between sites, we developed a containerized pipeline - all required software is included in this container and it runs independent of hardware. 
We adopted two container softwares
1. Docker – for all systems on computer with internet access
#### * Iremos utilizar Docker

2. Singularity – for unix only and for systems without internet-access

The software in the enigma-cnv container is:
•        PennCNV v1.0.5, used for CNV calling
•        R 3.3.1 + relevant R packages including iPsychCNV for visualization

### **1. Docker software**
You first need to install Docker:
https://docs.docker.com/get-docker/
Please follow the instructions to download docker for your particular platform. 
After download of docker software, please download the enigma-cnv:latest container by writing in the terminal:
```
docker pull bayramalex/enigma-cnv
```
Alternatively, find it on the ENIGMA-CNV github.

### **b. Scripts for running analysis** 
Please download theses scripts from github: 
https://github.com/ENIGMA-git/ENIGMA-CNV/CNVCalling/scripts/
ENIGMA_CNV_CNVProtocol_v2_singularity.sh or 
ENIGMA_CNV_CNVProtocol_v2_docker.sh [choose either singularity or docker- script]
ENIGMA-CNV_visualize_v1.R
compile_pfb_new.pl

### **c. Genotyping-chip-dependent files**
As mentioned,  genetic files and data are specific to genome-version. First wave of ENIGMA-CNV was based on the hg18 genome build. Most genotyping chips are now released with hg38 coordinates. In order of preference, we now aim for: hg38, hg19 and hg18. 
In addition to the correct genome version, appropriate files according to your genotyping chip must be used:
-The PFB-file (population frequency of B-allele file) - supplies the PFB information for each marker, and can give the chromosome coordinate information to PennCNV for CNV calling.
-The GCMODEL file - specifies the GC content of the 1Mb genomic region surrounding each marker (500kb each side) (calculated using the UCSC GC annotation file). ## Não temos
-the HMM-file -  tells the program what would be the expected signal intensity values for different copy number state, and what is the expected transition probability for different copy number states (For more information, see here: 
http://penncnv.openbioinformatics.org/en/latest/user-guide/input/)

### **Generating your own PFB-file**
For cohorts with more than 300 individuals with good quality genotyping data (for a suggestion good data, see footnote), we advise that you create your own PFB-file and GC-model files  – *this procedure is part of the script*.

### **The HMM-file**
There are several HMM-files:
hhall.hmm – can be used for all Illumina arrays (path, container: /opt/PennCNV- 1.0.5/lib/) ## Colocar esse caminho no pipeline

### **CNV of interest file based on genotyping chip**
We visualize a selected set of CNVs  - both with set and non-set boundaries. 
For this we need one of the following files:
CNVsofInterest_ENIGMA-CNV_hg18.csv
CNVsofInterest_ENIGMA-CNV_hg19.csv ## Nosso RG
CNVsofInterest_ENIGMA-CNV_hg38.csv
This list of CNVs is compiled from the paper from UK biobank: Kendall et al., 2017 (doi: 10.1016/j.biopsych.2016.08.014)  with a few additions.
Present at: https://github.com/ENIGMA-git/ENIGMA- CNV/CNVCalling/CNVsofInterests/

*Download according to your genome version and put in your Analysis-folder.*

### **Filtering genomic regions files**
Several genomic regions are known to harbor spurious CNV calls that might represent cell-line artifacts. 
Therefore, we remove CNVs in certain genomic regions overlapping >50%: 
-centromeric regions
-telomeric regions
-segmental duplication regions
-immunoglobulin regions 
You can find the corresponding files here: https://github.com/ENIGMA- git/ENIGMA-CNV/CNVCalling/filtergenomeregions/, e.g.:

*Download according to your genome version and put in your Analysis-folder.*

### **3. Call and visualize CNVs**

######## LER ANTES DE RODAR OFICIALMENTE

### **a. Adjust the script**
a. Please rename the script ENIGMA_CNV_CNVProtocol_v2.sh to ${Dataset}_ENIGMA_CNV_CNVProtocol_v2.sh [i.e. the name of your file with your dataset name [as stated below] as prefix]
b. edit the script by changing to your dataset-specific files as directed in section 0a [labeled “USER-INPUT needed] in the script.
d. Regarding deidentification – the protocol has an option to add an additional step of deidentification to the IDs – this leaves a key in the hands of the cohort submitting data. Adding this can give additional complications (if trying to track back to data). Please only tick this if required by your IRB-approvals. 
b. Run the analysis
Run the ${Dataset}_ENIGMA_CNV_CNVProtocol_v2.sh script in the terminal: bash ./${Dataset}_ENIGMA_CNV_CNVProtocol_v2.sh

After the run is done, please check the checklist found in ${Dataset}_visualize/${Dataset}_checklist.txt to see if things seem to make sense.
The CNV calling data is now ready for transfer (see E. below). 

## **Manipulação dos inputs**

- Identificamos os IDs com neuroimagem, mudamos para o ID do array no R (KeepFile.txt) e também fizemos o SexFile
	-Perda de 12 probandos que não tinham o ID da sonda no espelho
- Pegamos os VCF individuais nomeados como INPD_*_subjID com as colunas: ID	CHROM 	POS	subjID.B Allele Freq	subjID.Log R Ratio
	- Modificamos o nome com o seguinte comando para ficar só o subjID da sonda:
```
ls INPD* | xargs -n 1 | sed 's/^[^.]*.//' > new_filename.txt
files=(*)
i=0
while read -r new_name; do
  mv "${files[$i]}" "$new_name"
  (( i++ ))
done < new_filename.txt
```
Esses arquivos ficaram no Analysis_Folder/data_neuroimaging

- Selecionamos do julia_data (todos os individuos) os arquivos que estavam no KeepFile.txt, que são os que realmente serão usados com select_ids_neuroimaging.sh:

```
#!/bin/bash

# Set the source and destination directories
source_dir="/mnt/genetica_1/JuliaArendt/ENIGMA-CNV_Analysis/ENIGMA-CNV/CNVCalling/Analysis_Folder/julia_data"
destination_dir="/mnt/genetica_1/JuliaArendt/ENIGMA-CNV_Analysis/ENIGMA-CNV/CNVCalling/Analysis_Folder/data_neuroimaging"

# Name list file
name_list_file="../KeepFile.txt"

# Check if the name list file exists
if [ ! -e "$name_list_file" ]; then
    echo "Name list file not found: $name_list_file"
    exit 1
fi

# Loop through each line in the name list file and move files
while IFS= read -r filename; do
    # Trim leading and trailing whitespaces from the filename
    filename=$(echo "$filename" | tr -d '[:space:]')

    # Check if the file exists in the source directory
    if [ -e "$source_dir/$filename" ]; then
        # Move the file to the destination directory
        mv "$source_dir/$filename" "$destination_dir/$filename"
        echo "Moved: $filename"
    else
        echo "File not found: $filename"
    fi
done < "$name_list_file"
```
	- Perda de mais 20 probandos não sabemos por que (total = 708)

- Modificamos as colunas com manipulate_columns.sh no data_neuroimaging:

```
#!/bin/bash

# Set the directory containing your files
files_dir="./"

# List of words to be substituted
word1="ID"
word2="CHROM"
word3="POS"
word4="B Allele Freq" 
word5="Log R Ratio"

# List of corresponding new words
new_word1="Name"
new_word2="Chromosome"
new_word3="Position"
new_word4="BAF" # Mudamos essas 2 últimas para tirar o espaço na palavra para quando fossemos inverter as colunas (como está nas instruções) não desse trabalho
new_word5="LRR"

# Loop through all files in the specified directory
for file in "$files_dir"/*; do
    # Check if the file is a regular file
    if [ -f "$file" ]; then
        # Substitute the words using sed
        sed -i "s/$word1/$new_word1/g; s/$word2/$new_word2/g; s/$word3/$new_word3/g; s/$word4/$new_word4/g; s/$word5/$new_word5/g" "$file"
        echo "Substituted words in: $file"
    fi
done
```
- Em seguida, copiei um ID aleatório para o Analysis folder para excluir as 2 últimas colunas (BAF e LRR) e fazer o SNP-position-file (Name Chromosome Position) com o seguinte comando:

```
awk -F'\t' '{print $1"\t"$2"\t"$3}' 203273230127_R08C02 > SNP-position-file # Importante manter o tab depois de cada coluna (SNP-position-file = 688033 linhas/variantes)
```

- Agora precisamos tirar as colunas Chromosome e Position e inverter BAF e LRR, mas em todos os arquivos em data_neuroimaging

```
#!/bin/bash

# Specify the input directory and output directory
input_directory="./"
output_directory="./neuroimaging_data_LRR_BAF"

# Create the output directory if it doesn't exist
mkdir -p "$output_directory"

# Loop through all files in the input directory
for file in "$input_directory"/2*; do
    # Extract the filename without the path
    filename=$(basename "$file")

    # Define the output file path
    output_file="$output_directory/$filename"

    # Use awk to extract and reorder columns
    awk -F'\t' '{print $1"\t"$5"\t"$4}' "$file" > "$output_file"
done
```
- Para mudar de LRR e BAF para o nome extenso: 

```
#!/bin/bash

# Set the directory containing your files
files_dir="./"

# List of words to be substituted
word1="LRR"
word2="BAF"

# List of corresponding new words
new_word1="Log R Ratio"
new_word2="B Allele Freq"

# Loop through all files in the specified directory
for file in "$files_dir"/*; do
    # Check if the file is a regular file
    if [ -f "$file" ]; then
        # Substitute the words using sed
        sed -i "s/$word1/$new_word1/g; s/$word2/$new_word2/g; s/$word3/$new_word3/g; s/$word4/$new_word4/g; s/$word5/$new_word5/g" "$file"
        echo "Substituted words in: $file"
    fi
done
```

- Como 32 individuos foram excluídos, atualizar o KeepFile com ls 2* | xargs -n 1 > novo_Keep.txt dentro da pasta neuroimaging_data_LRR_BAF e movê-lo para Analysis_Folder 


## **Rodando a pipeline de CNVS**


```
#!/bin/bash

set -xv #print commands
exec > >(tee -i ENIGMA-CNV_logfile.txt) # write output to logfile # USER-input - if running on a windows machine, this may need to be commented out
exec 2>&1 # ensures that all std/stderr are also written to log-file # USER-input - if running on a windows machine, this may need to be commented out

#set -xv #print commands

# date and start the clock
date=$(date +%y%m%d)
START=$(date +%s)

#: Document: ENIGMA-CNV_CNVProtocol_v2_docker.sh
#: Title: ENIGMA CNV Working group – standardized protocol for calling and filtering CNVs
#: Purpose:
	# Call raw CNVs on Illumina samples and Affymetrix samples [the latter needs preprocessing first]
	# Filtering of raw CNVs
	# Visualization of selected CNVs of Interest and CNVs>50kb
# Date, v1: 2015-07-06
# Date, v2: 2022-02
# Authors, v1: Ida Elken Sønderby, Omar Gustafsson. Input and tested by Allan McRae and Nicola Armstrong
# Versions: v1.1. Minor bugs # v1.2: Changed the protocol so that it can be run as a bash-script. Implemented Rscripts for downstream filtering of spurious CNVs # v1.3: Changed incorporation of keep and remove-file
# Authors, v2: Ida Elken Sønderby.
# v2.0: Altered to use software containers, included visualization in script

################################
#### INFO prior to starting ####
################################

# OBS – only data from a single chip can be run - if data is available from more than one type of chip, please run the protocol separately for each different chip set

# 1. Please read and follow the instructions in the word-document "Instructions, ENIGMA-CNV working group, v2.0" before commencing analysis, which includes:
# - identifying appropriate files for CNV calling
# - downloading and installing required software and scripts
# - adjust the names and directory-instructions in the script itself as directed in 0a below [labeled USER-INPUT needed].

# 2. Run the script:
    # a. Prior to this, you may need to change file permission in terminal: chmod +u+x ${Dataset}_ENIGMA_CNV_CNVProtocol_v2_docker.sh
    # b. Run script in terminal: bash ./${Dataset}_ENIGMA_CNV_CNVProtocol_v2_docker.sh

# 3. After the run, the printed log will be present in the Analysis-folder: ENIGMA-CNV_logfile.txt (in case something goes wrong, you can hopefully retrack it from here, please check for error-messages)

# 4. After, the run, please check the checklist found in ${Dataset}_visualize/${Dataset}_checklist.txt to see if things make sense.

# 5. Please, send the compressed files in the Analysis-dir: ${Dataset}_visualize.tar.gz-folder to: enigmacnvhelpdesk@gmail.com

# Please address any questions to: enigmacnvhelpdesk@gmail.com

##################################################
## 0a. USER-files and input - USER input needed ##
##################################################

# General info
declare Dataset="BHRCS" # USER-INPUT - Replace test with the name of your dataset
declare ResponsiblePI="Santoro" # USER-INPUT - Replace OleAndreassen with the name of the PI [no spaces, please]
declare Email_PI="santorogen@gmail.com" # USER-INPUT - Replace o.a.andreassen@medisin.uio.no with the e-mail of your PI to enable us to contact the PI directly in case of questions.
declare Analyst="Arendt" # USER-INPUT - Replace IdaSoenderby with the name of the analyst [no spaces, please]
declare Email_Analyst="julia.arendt@unifesp.br" # USER-INPUT - Replace e-mail with the e-mail of the analyst to enable us to contact the analyst directly in case of questions.
declare ANALYSISDIR="/mnt/genetica_1/JuliaArendt/ENIGMA-CNV_Analysis" # USER-INPUT - replace with the full path for the Analysis-folder on your computer/server
declare deidentify="NA" # USER-INPUT - only change to YES if deidentification is really important to you, otherwise keep NA

# a. Genetic information files - input specific to your genotyping chip
declare Chip="Illumina_Global_Screening_Array" # USER-INPUT - genotyping chip: Note - if your dataset is Affymetrix, make sure to put in "Affy6" - otherwise the script will not run appropriately
declare Chipversion="V3" # USER-INPUT - version of chip (if applicable) [set as NA if non-applicable]
#declare ILLUMINAREPORTDIR=""  # USER-INPUT - absolute path to your Illumina-Final Report-file. Can be replaced with ${ANALYSISDIR} if you placed it there. If you already have LRRBAF-files, you do not need to fill it out
declare IlluminaReport="" # USER-INPUT - name of your Illumina-Final Report-file ## Apaguei o que tinha aqui e em cima pq temos LRR/BAF
declare SNPPosFile="SNP-position-file" # USER-INPUT
	# The snp-position-file is a tab-delimited file with the positions of the SNPs on your chip, containing at least the columns below (with the exact (!) headers). These columns are for instance present in the SNP-map-file [.map] generated together with your IlluminaReport or the Illumina manifest file [.bpm]:
	# Name	Chromosome	Position
	# rs1000000	12	126890980
	# rs1000002	3	183635768
	# rs10000023	4	95733906
declare genomeversion="hg19" # USER-INPUT - Genome-version of SNPpositionfile. It is important that this is correct - otherwise CNV calls and visualization become wrong.
# IF BAF-LRR-files already present for your dataset (for instance Affy-users)
declare BAFLRRpresent="yes" # USER-INPUT - only put "yes" if you have BAF-LRR-files ready available - otherwise, leave the "no"
if [ $BAFLRRpresent = "yes" ]
    then
    declare LRRBAFDIR="${ANALYSISDIR}/BAF_LRR/" # USER-INPUT -  directory with PennCNV-inputfiles
# NOTE - if the lrr-baf-folder has several subfolders with lrrbaf-files, the script may run into issues... If so, please contact the helpdesk
    else
    mkdir ${ANALYSISDIR}/LRRBAF_Files/
declare LRRBAFDIR="${ANALYSISDIR}/LRRBAF_Files/" # Note - this will place the PennCNV input files in a subfolder of  your ENIGMA-CNV analysis folder. These take up quite a lot of space - if you wish them to be placed elsewhere, write the full path of the wanted folder (USER-INPUT)
fi
# Post/prefixes in LRR-BAF-files
declare postscript="" # e.g ".penn" # USER-INPUT - IF converting from IlluminaFinalReport to LRRBAF-files with this script, this should be left empty. LRR-BAF files from previous convertions may have a postfix (e.g. ".penn", "_lrrbaf"), please input this here (USER-INPUT).
#declare Affyprescript="" # e.g. "gw." # USER-INPUT - Affymetrix files names often get this prescript, which needs to be removed to couple samples. Replace with different prescript if relevant (USER-INPUT).

# b. Cohort-generated files
declare SexFile="${ANALYSISDIR}/SexFile" # USER-INPUT - absolute path to your sex-file
# If sex for a sample is not provided in sexfile, or if --sexfile is not specified, PennCNV will try to predict the gender of the sample. It is highly recommended to provide a sexfile [saves time].
declare gendermissing="0" # USER-INPUT, number of individuals with gender missing in sexfile
#declare RelativeFile="${ANALYSISDIR}/DupsRelatives.txt" # USER-INPUT - absolute path to your relative-file
declare KeepFile="" # USER-INPUT # absolute path to your KeepFile.txt - leave empty if you do not have a KeepFile
#declare RemoveFile="${ANALYSISDIR}/RemoveFile.txt" # USER-INPUT # absolute path to your RemoveFile.txt - leave empty if you do not have a RemoveFile

### c. Genotyping-chip dependent information -
declare HMMname="/opt/PennCNV-1.0.5/lib/hhall.hmm" # USER-INPUT - Please replace with the correct HMM-file. examples:  /opt/PennCNV-1.0.5/lib/hhall.hmm; for /opt/PennCNV-1.0.5/libhh550.hmm /opt/PennCNV-1.0.5/affy/libgw6/affygw6.hmm

# IF your dataset has more than 300 individuals, you can generate your own PFB-file based on the frequency in your dataset
#declare NoofIndividuals="" # USER-INPUT - e.g "1000" - No of individuals (e.g. "300") to be used for generating PFB-file and GCC-file (must be at least 300 individuals of good quality). Leave empty if you want to use all your individuals. The more individuals you use, the more precise the estimate becomes but the longer it will take to generate the PFB-model. For the NORMENT dataset, generating a PFB- and GCMODEL-file for the OmniExpress12v1.0 containing 730,525 markers using 1000 individuals took ~90 min on a Mac laptop with a 2.53 GHz Intel Core i5 processor and 4 GB of working memory

# IF your dataset contains less than 300 individuals, you need to use a generic version of the PFB-file. Please confer with the ENIGMA-CNV working group and put in the correct names of the files. # NOTE - IF more than 300 individuals, these files will be generated later and named "${Dataset}_${genomeversion}.pfb" (likewise for GCName)
#declare PFB="hhall.${genomeversion}" # ONLY USER-INPUT for those with <300 individuals -  replace with the correct PFB & GCMODEL (note without extensions)
declare PFBname="BHRCS_hg19.pfb"  # no input
declare GCname="BHRCS_hg19.gcmodel" # no input

#####################################################################
## 0b. State cut-offs and predefined files/folders - NO user input ##
#####################################################################

# a. List of input and output-files from ENIGMA CNV calling protocol
declare List_preQC="${Dataset}_ListofInputFiles_preQC.txt"
declare List_postQC="${Dataset}_ListofInputFiles_postQC.txt"

# b. Cut-offs, autosomal chromosomes filtering
declare LRR_SD=0.40
declare BAF_drift=0.02
declare WF=0.05
declare NoofSNPs=15
declare MergeFraction=0.30 # depending on chip, this may need adjustment - this was appropriate for Illumina OmniExpress where 0.2 was too low.
declare MinQueryFrac=0.5 # define the overlap with the regions necessary to be excluded

# c. Cut-offs, X chromosome filtering
# For the X-chromoxome, only small CNVs are removed and CNVs merged whereas BAF-drift and WF, LRR_SD are skipped (by setting them abnormally high) to not filter based on X-chromosome only.
declare NoofSNPs_X=15
declare MergeFraction_X=0.30
declare LRR_SD_X=0.99
declare BAF_drift_X=0.99
declare WF_X=0.99

### d. Parameters needed for visualization
declare VISUALIZEDIR=${ANALYSISDIR}/${Dataset}_visualize/
mkdir ${VISUALIZEDIR}
declare CNVofInterestFile="CNVsofInterest_ENIGMA-CNV_hg19.csv" # appropriate CNVsofInterest-file for visualization (changes according to genome version)
Overlapref=0.3 # How large a proportion of the CNV is overlapping with the CNVofInterest?
OverlapMin=0.35 # Minimum overlap
OverlapMax=5 # Maximum overlap

### e. Output directories (to not clutter Analysis-dir)
declare OUTDIR=${ANALYSISDIR}/${Dataset}_output
mkdir ${OUTDIR}

##################
#### ANALYSIS ####
##################

###########################
## The protocol in short ##
###########################

#### A. Generate input-files for PennCNV
	## Step 1: Generate inputfiles from Illumina FinalReport files
	## Step 2: Make a list of sample inputfiles for PennCNV

#### B. Select helper files
	## Step 1. Select helper-files

#### C: Do CNV Calling
	## Step 1: Call CNVs on autosomal chromosomes with GC-adjustment
	## Step 2: Call CNVs on X-chromosome with GC-adjustment

#### D: First round of QC
	## Step 1: QC on CNVs, 1st round. Obtain summary statistics
	## Step 2. Merge CNVs
	## Step 3. Remove spurious CNV calls in specific genomic regions
	## Step 4. Obtain summary statistics for the QCed dataset
	## Step 5. Removing QC'ed individuals from sumout-lists

#### E: Deidentification

#### F: CNV visualization

#### G: Data transfer
	## Step 1: Make folder for file-transfer

#########################################
## A. Generate input-files for PennCNV ##
#########################################

###############################################################
# Step 1: Generate inputfiles from Illumina FinalReport files #
###############################################################

# These steps convert the Illumina Report files into separate files for each individual with intensity data for each individual sample
if [ $BAFLRRpresent != "yes" ]
then
# a. Split IlluminaReport into one file for each individual sample
sudo docker run -v ${ILLUMINAREPORTDIR}:/illuminadir -v ${LRRBAFDIR}:/lrrbafdir  bayramalex/enigma-cnv:latest /opt/PennCNV-1.0.5//split_illumina_report.pl --prefix /lrrbafdir/ /illuminadir/${IlluminaReport}
fi
	#  "-prefix"	used to specify the prefix of the file name (in this case, save the file to the ${ANALYSISDIR}/LRRBAF_Files/-directory).
	#  "-comma"	can be added to the command in case the genotyping center provide the report file in comma-delimited format rather than tab-delimited format

	# Output of this file is a file for each individual named after the “Sample ID” field (e.g. NORMENT1) in the original IlluminaReport file with the following (tab-delimited) format:
		# Name	NORMENT1.Log R Ratio	NORMENT1.B Allele Freq
		# SNP1	-0.0038	0.0161
		# SNP2	0.0073	0.9943
		# rs10000023	-0.0307	0.0026
		# 	etc…

##########################################################
## Step 2: Make a list of sample inputfiles for PennCNV ##
##########################################################

# a. Generate a list-file with one individual on each line with the full path for its inputfile (i.e. the output from step1)
find ${LRRBAFDIR} -type f -name '*' > ${OUTDIR}/${List_preQC}

# b. Keep individuals based on your KeepFile.txt
if [ -f "${KeepFile}" ] # If Keepfile exist
   then # keep only those samples
		awk '{print $1,$1}' ${OUTDIR}/${List_preQC} | awk -v postscript=${postscript} -v prescript=${prescript} '{gsub(/.+\//,"\/",$1); gsub(/\//,"",$1); gsub(postscript,"",$1); gsub(prescript,"",$1); print}' | awk 'FNR==NR {k[$1];next} {if ($1 in k) {print $2} else {}}' ${KeepFile} -  >${OUTDIR}/${List_postQC}
   else # keep everyone
   cp ${OUTDIR}/${List_preQC} ${OUTDIR}/${List_postQC}
fi

# c. -and remove individuals based on your removefile
if [ -f "$RemoveFile" ] # if Removefile exist
   then # Remove samples
        	awk '{print $1,$1}' ${OUTDIR}/${List_postQC} | awk -v postscript=${postscript} -v prescript=${prescript} '{gsub(/.+\//,"\/",$1); gsub(/\//,"",$1); gsub(postscript,"",$1); gsub(prescript,"",$1); print}' | awk 'FNR==NR {k[$1];next} {if ($1 in k) {} else {print $2}}' ${RemoveFile} - >tmp
	mv tmp ${OUTDIR}/${List_postQC}
fi

## Change input-file list so it is compatible with containers

# Adjust to make input for container
cp ${OUTDIR}/${List_postQC} ${OUTDIR}/${List_postQC}_adj
sed -i -e 's/.*\//\/lrrbafdir\//' ${OUTDIR}/${List_postQC}_adj
# NOTE - if files come from several folders for lrr-baf-files, this command may run into issues... Please contact the helpdesk

#######################################
# Step 3: Generate PFB & GCMODEL-file #
#######################################

	# Description: The PFB-file (population frequency of B-allele file) supplies the PFB information for each marker, and gives the chromosome coordinates information to PennCNV for CNV calling. It is a tab-delimited text file with four columns, representing marker name, chromosome, position and PFB values.
	# The script compile_pfb_new.pl compiles a PFB file from multiple signal intensity files containing BAF values
	# The script cal_gc_snp.pl calculates GC content surrounding each marker within specified sliding window, using the UCSC GC annotation file.

#declare length_fileinput=`wc -l ${OUTDIR}/${List_postQC} | awk '{print $1}'` # number of individuals in

#if [ ${length_fileinput} -gt 299 ]
#	then
#  echo "more than 299 individuals, a PFB-file for your dataset will be created"
#	declare PFBname="${Dataset}_${genomeversion}.pfb" # PFB-file
#	declare GCname="${Dataset}_${genomeversion}.gcmodel" # GC-model-file
#	declare List_Helperfiles="${Dataset}_ListofInputFiles_PFB.txt" # files used for input for PFB

#	if [ -n "${NoofIndividuals}" ] ## variable exists and is not empty
#	then
#		echo "number of individuals used for generating PFB-file: ${NoofIndividuals}"
#	else
#		declare NoofIndividuals=${length_fileinput}
#		echo "number of individuals used for generating PFB-file: ${NoofIndividuals}"
#	fi
		# a. Generate a list of x individuals randomly selected from your file
#		export NoofIndividuals; awk 'BEGIN {p=ENVIRON["NoofIndividuals"]; srand();} {a[NR]=$0} END{for(i=1; i<=p; i++){x=int(rand()*NR) + 1; print a[x]}}' ${OUTDIR}//${List_postQC} >${OUTDIR}/${List_Helperfiles}
		# b. Generation of PFB-file
#		perl compile_pfb_new.pl --listfile ${OUTDIR}/${List_Helperfiles} -snpposfile ${ANALYSISDIR}/${SNPPosFile} --output ${OUTDIR}/${PFBname}

			# --output <file>             specify output file (default: STDOUT)
			# --snpposfile <file>         a file that contains Chr and Position for SNPs
			# --listfile <file>           a listfile that contains signal file names

	#cp BHRCS_output/BHRCS_hg19.pfb ${ANALYSISDIR}

		
		# c. Generation of GC-model file
#		sudo docker run -v ${ANALYSISDIR}:/analysisdir bayramalex/enigma-cnv:latest cp /opt/PennCNV-1.0.5/gc_file/${genomeversion}.gc5Base.txt.gz /analysisdir/ # copy gcbase-file over from penncnv-dir to be able to unzip
#		gunzip ${ANALYSISDIR}/${genomeversion}.gc5Base.txt.gz
#		sudo docker run -v ${OUTDIR}:/outdir -v ${ANALYSISDIR}:/analysisdir bayramalex/enigma-cnv:latest /opt/PennCNV-1.0.5/cal_gc_snp.pl /analysisdir/${genomeversion}.gc5Base.txt /analysisdir/${PFBname} -output /analysisdir/${GCname}
#		rm ${ANALYSISDIR}/${genomeversion}.gc5Base.txt
#else
#	echo "less than 300 individuals, please choose a generic PFB-file after conferring with the ENIGMA-CNV working group"
#fi

#######################
## C: Do CNV Calling ##
#######################

# Autosomal and X-chromosome CNVs are called separately

## Step 1: Call CNVs on autosomal chromosomes with GC-adjustment
sudo docker run -v ${OUTDIR}:/outdir -v ${ANALYSISDIR}:/analysisdir -v ${LRRBAFDIR}:/lrrbafdir bayramalex/enigma-cnv:latest /opt/PennCNV-1.0.5/detect_cnv.pl -test -hmm ${HMMname} -pfb /analysisdir/${PFBname} -gcmodel /analysisdir/${GCname} -list outdir/${List_postQC}_adj --confidence -out outdir/${Dataset}.auto.raw  -log outdir/${Dataset}.auto.raw.log

	# --test		tells the program to generate CNV calls
	# --confidence	calculate confidence for each CNV
	# --hmmfile and –pfbfile	provides the files for hmm- & pfb-files
	# --gcmodel	implements a wave adjustment procedure for genomic waves which reduce false positive calls. It is recommended when if a large fraction of your samples have waviness factor (WF value) less than -0.04 or higher than 0.04
	# -list		provides the file with the list of samples you want called
	# –log		tells the program to write log information
	# --out		raw cnvcalls for all input individuals

# Note: This took ~48 hours to run for the MAS sample which comprises 925 individuals and 1,8 mio marker on the Affymetrix genomewide array 6.0-platform.


## Step 2: Call CNVs on X-chromosome with GC-adjustment

# a. Make sex-file for your inputfile
awk '{print $1,$1}' ${OUTDIR}${List_postQC}_adj | awk -v postscript=${postscript} -v prescript=${prescript} '{sub(/.+\//,"\/",$1); sub(/\//,"",$1); gsub(postscript,"",$1); gsub(prescript,"",$1); print}' | awk 'FNR==NR {k[$1]=$2;next} {if ($1 in k) {print k[$1] "\t" $2}}' - ${SexFile} >${OUTDIR}/${Dataset}_SexFile.txt

    # If sex for a sample is not provided in sexfile, or if --sexfile is not specified, PennCNV will try to predict the gender of the sample. It is highly recommended to provide a sexfile [saves time].
    # This command couples the full path of each individual input-file to sex with the following end-format:
    # 	/Volumes/CNV_Calling/Analysis/PennCNV_InputFiles/TOP1	female
    # 	/Volumes/CNV_Calling/Analysis/PennCNV_InputFiles/TOP2	male
    # 	etc

# b. Call CNVs on X-chromosome
sudo docker run -v ${OUTDIR}:/outdir -v ${ANALYSISDIR}:/analysisdir -v ${LRRBAFDIR}:/lrrbafdir bayramalex/enigma-cnv:latest /opt/PennCNV-1.0.5/detect_cnv.pl -test -hmm ${HMMname} -pfb /analysisdir/${PFBname} -gcmodel /analysisdir/${GCname} -list outdir/${List_postQC}_adj --confidence -out outdir/${Dataset}.X.raw  -log outdir/${Dataset}.X.raw.log --chrx --sexfile outdir/${Dataset}_SexFile.txt
	# --chrx		specifies that x-chromosome should be called
	# --sexfile	provides the sexfile data for the calling

# Note: This took ~48 hours to run for the MAS sample which comprises 925 individuals and 1,8 mio marker on the Affymetrix genomewide array 6.0-platform on a ??? specify computer etc...

######################
## D. Downstream QC ##
######################

###########
# OUTLINE #
###########

# The raw CNV calls need to have calls from low quality samples eliminated and low quality calls need to be removed.

# AUTOSOMAL CNVs:
# Step 1: QC on CNVs, 1st round. Obtain summary statistics
# Step 2. Merge CNVs
# Step 3. Removing spurious CNV calls in specific genomic regions
# Step 4. Obtain summary statistics for the QCed dataset

# X-chromosomal CNVs:
# Step 1: QC on CNVs, 1st round. Obtain summary statistics
# Step 2. Merge CNVs
# Step 3. Remove spurious CNV calls in specific genomic regions
# Step 4. Obtain summary statistics for the QCed dataset

# Common:
# Step 5. Removing QC'ed individuals from sumout-lists
# Step 6: Putting together checklist-file


####################
## AUTOSOMAL CNVs ##
####################

############################################################
# Step 1: QC on CNVs, 1st round. Obtain summary statistics #
############################################################

# The filter_cnv.pl program identifies low-quality samples from a genotyping experiment and eliminates them from future analysis. This analysis requires the output LOG file from CNV calling in addition to the raw cnv-file.

# a. obtain summary statistics for dataset
sudo docker run -v ${OUTDIR}:/outdir -v ${ANALYSISDIR}:/analysisdir   bayramalex/enigma-cnv:latest /opt/PennCNV-1.0.5/filter_cnv.pl outdir/${Dataset}.auto.raw --qclrrsd ${LRR_SD} --qcbafdrift ${BAF_drift} --qcwf ${WF} --numsnp ${NoofSNPs} --qclogfile outdir/${Dataset}.auto.raw.log --qcpassout outdir/${Dataset}.auto.passout --qcsumout outdir/${Dataset}.auto.sumout --out outdir/${Dataset}.auto.flr

echo "Finished first filtering of autosomal CNVs"

######################
# Step 2: Merge CNVs #
######################


# i. 1st time
sudo docker run -v ${OUTDIR}:/outdir -v ${ANALYSISDIR}:/analysisdir  bayramalex/enigma-cnv:latest /opt/PennCNV-1.0.5/clean_cnv.pl combineseg --fraction ${MergeFraction} --bp --signalfile /analysisdir/${PFBname} outdir/${Dataset}.auto.flr --output outdir/${Dataset}.autosomal.flr_mrg1

# ii. This command ensures that CNVs are getting merged until there are no more CNVs to merge within the defined distance
{
i=1
while [ ${i} -lt 999 ]; do
declare j=$(($i+1));
sudo docker run -v ${OUTDIR}:/outdir -v ${ANALYSISDIR}:/analysisdir bayramalex/enigma-cnv:latest /opt/PennCNV-1.0.5//clean_cnv.pl combineseg --fraction ${MergeFraction} --bp --signalfile /analysisdir/${PFBname} outdir/${Dataset}.autosomal.flr_mrg${i} --output outdir/${Dataset}.autosomal.flr_mrg${j}

declare length1=`awk 'END {print NR}' ${OUTDIR}/${Dataset}.autosomal.flr_mrg${i}`
declare length2=`awk 'END {print NR}' ${OUTDIR}/${Dataset}.autosomal.flr_mrg${j}`
if [ ${length1} -eq ${length2} ]
then
cp ${OUTDIR}/${Dataset}.autosomal.flr_mrg${j} ${OUTDIR}/${Dataset}.auto.flr_mrg_final
break;
fi
i=${j};
done
}
echo "Finished merging of autosomal CNVs"

# remove the resulting merging files
if test -f "${OUTDIR}/${Dataset}.autosomal.flr_mrg[0-9]"; then
	rm ${OUTDIR}/${Dataset}.auto.flr_mrg[0-9]
fi
if test -f "${OUTDIR}/${Dataset}.autosomal.flr_mrg[0-9][0-9]"; then
	rm ${OUTDIR}/${Dataset}.auto.flr_mrg[0-9][0-9]
fi

#################################################################
# Step 3. Remove spurious CNV calls in specific genomic regions #
#################################################################
# Several genomic regions are known to harbor spurious CNV calls that might represent cell-line artifacts. Therefore, we remove CNVs in certain genomic regions:

# a. Identify overlapping CNVs

# i. Identify CNVs with overlap to centromeric, telomeric, segmentalduplication and immunoglobulin regions
for i in centro telo segmentaldups immuno;
do
	sudo docker run -v ${OUTDIR}:/outdir -v  ${ANALYSISDIR}:/analysisdir bayramalex/enigma-cnv:latest /opt/PennCNV-1.0.5/scan_region.pl outdir/${Dataset}.auto.flr_mrg_final analysisdir/${i}_${genomeversion}.txt -minqueryfrac ${MinQueryFrac} >${OUTDIR}/${Dataset}.auto.${i};
echo "${i} is done";
done

b. Remove CNVs with overlap to spurious regions

# i. Make list
for i in centro telo segmentaldups immuno;
do
	cat ${OUTDIR}/${Dataset}.auto.${i};
done >${OUTDIR}/${Dataset}_SpuriousCNVs_exclude

# ii. Remove CNVs in spurious  regions
grep -Fv -f ${OUTDIR}/${Dataset}_SpuriousCNVs_exclude ${OUTDIR}/${Dataset}.auto.flr_mrg_final >${OUTDIR}/${Dataset}.auto.flr_mrg_spur
echo "Finished removal of spurious regions for autosomal CNVs"

##########################################################
# Step 4. Obtain summary statistics for the QCed dataset #
##########################################################

# a. obtain summary statistics for dataset
sudo docker run -v ${OUTDIR}:/outdir -v  ${ANALYSISDIR}:/analysisdir bayramalex/enigma-cnv:latest /opt/PennCNV-1.0.5/filter_cnv.pl outdir/${Dataset}.auto.flr_mrg_spur -qclrrsd ${LRR_SD} --qcbafdrift ${BAF_drift} --qcwf ${WF} -numsnp ${NoofSNPs} --qclogfile outdir/${Dataset}.auto.raw.log --qcpassout outdir/${Dataset}.auto.passout_QC --qcsumout outdir/${Dataset}.auto.sumout_QC --out outdir/${Dataset}.auto.flr_QC

#######################
# X-chromosomal CNVs ##
#######################

###########################################################
# Step 1: QC on CNVs, 1st round. Obtain summary statistics #
###########################################################

# a. obtain summary statistics for dataset
sudo docker run -v ${OUTDIR}:/outdir -v ${ANALYSISDIR}:/analysisdir   bayramalex/enigma-cnv:latest /opt/PennCNV-1.0.5/filter_cnv.pl outdir/${Dataset}.X.raw --qclrrsd ${LRR_SD_X} --qcbafdrift ${BAF_drift_X} --qcwf ${WF_X} --numsnp ${NoofSNPs_X} --qclogfile outdir/${Dataset}.X.raw.log --qcpassout outdir/${Dataset}.X.passout --qcsumout outdir/${Dataset}.X.sumout --out outdir/${Dataset}.X.flr --chrx

echo "Finished first filtering of X-chromosomal CNVs"

######################
# Step 2. Merge CNVs #
######################

# i. 1st time
sudo docker run -v ${OUTDIR}:/outdir -v ${ANALYSISDIR}:/analysisdir  bayramalex/enigma-cnv:latest /opt/PennCNV-1.0.5/clean_cnv.pl combineseg --fraction ${MergeFraction_X} --bp --signalfile /analysisdir/${PFBname} outdir/${Dataset}.X.flr --output outdir/${Dataset}.X.flr_mrg1

# ii. This command ensures that CNVs are getting merged until there are no more CNVs to merge within the defined distance
{
i=1
while [ ${i} -lt 999 ]; do
declare j=$(($i+1));
sudo docker run -v ${OUTDIR}:/outdir -v ${ANALYSISDIR}:/analysisdir bayramalex/enigma-cnv:latest /opt/PennCNV-1.0.5//clean_cnv.pl combineseg --fraction ${MergeFraction_X} --bp --signalfile /analysisdir/${PFBname} outdir/${Dataset}.X.flr_mrg${i} --output outdir/${Dataset}.X.flr_mrg${j}

declare length1=`awk 'END {print NR}' ${OUTDIR}/${Dataset}.X.flr_mrg${i}`
declare length2=`awk 'END {print NR}' ${OUTDIR}/${Dataset}.X.flr_mrg${j}`
if [ ${length1} -eq ${length2} ]
then
cp ${OUTDIR}/${Dataset}.X.flr_mrg${j} ${OUTDIR}/${Dataset}.X.flr_mrg_final
break;
fi
i=${j};
done
}
echo "Finished merging of X-chromosomal CNVs"

# remove the resulting merging files
if test -f "${OUTDIR}/${Dataset}.X.flr_mrg[0-9]"; then
	rm ${OUTDIR}/${Dataset}.X.flr_mrg[0-9]
fi
if test -f "${OUTDIR}/${Dataset}.X.flr_mrg[0-9][0-9]"; then
	rm ${OUTDIR}/${Dataset}.X.flr_mrg[0-9][0-9]
fi

###################################################################
# Step 3. Removing spurious CNV calls in specific genomic regions #
###################################################################

# Several genomic regions are known to harbor spurious CNV calls that might represent cell-line artifacts. Therefore, we remove CNVs in certain genomic regions:

#  a. Identify overlapping CNVs

# i. Identify CNVs with overlap to centromeric, telomeric and  segmentalduplication regions
for i in centro telo segmentaldups immuno;
do
sudo docker run -v ${OUTDIR}:/outdir -v  ${ANALYSISDIR}:/analysisdir bayramalex/enigma-cnv:latest /opt/PennCNV-1.0.5/scan_region.pl outdir/${Dataset}.X.flr_mrg_final /analysisdir/${i}_${genomeversion}.txt -minqueryfrac ${MinQueryFrac} >${OUTDIR}/${Dataset}.X.${i};
echo "${i} is done";
done

# b. Remove CNVs with overlap to spurious regions

# i. Make list
for i in centro telo segmentaldups immuno;
do cat ${OUTDIR}/${Dataset}.X.${i};
done >${OUTDIR}/${Dataset}.X_SpuriousCNVs_exclude

# ii. Remove CNVs in spurious regions
grep -Fv -f ${OUTDIR}/${Dataset}.X_SpuriousCNVs_exclude ${OUTDIR}/${Dataset}.X.flr_mrg_final >${OUTDIR}/${Dataset}.X.flr_mrg_spur
echo "Finished removal of spurious regions for autosomal CNVs"

##########################################################
# Step 4. Obtain summary statistics for the QCed dataset #
##########################################################

# a. Removing CNVs from the X-chr based on autosomal individual QC
# Some individuals were removed in the autosomal QC -  the CNVs related to these should be removed from the X-chromosomal dataset before getting the final summary statistics.
awk 'FNR==NR {a[$1]; next} {if ($5 in a) {print}}' ${OUTDIR}/${Dataset}.auto.passout_QC ${OUTDIR}/${Dataset}.X.flr_mrg_spur >${OUTDIR}/${Dataset}.X.flr_mrg_spur_onlypass

# b. obtain summary statistics for dataset
sudo docker run -v ${OUTDIR}:/outdir -v  ${ANALYSISDIR}:/analysisdir bayramalex/enigma-cnv:latest /opt/PennCNV-1.0.5/filter_cnv.pl outdir/${Dataset}.X.flr_mrg_spur_onlypass -qclrrsd ${LRR_SD_X} --qcbafdrift ${BAF_drift_X} --qcwf ${WF_X} -numsnp ${NoofSNPs_X} -qclogfile outdir/${Dataset}.X.raw.log -qcpassout outdir/${Dataset}.X.passout_QC -qcsumout outdir/${Dataset}.X.sumout_QC -out outdir/${Dataset}.X.flr_QC --chrx

##########################
## BOTH AUTOSOMAL AND X ##
##########################

########################################################
# Step 5. Removing QC'ed individuals from sumout-lists #
########################################################

# The sum-out list includes all individuals in the log-file. We wish to make a sum-out list only with individuals passing QC

# a. Removing individuals
for i in ${Dataset}.auto ${Dataset}.X;
do
	awk 'BEGIN {FS=OFS="\t"} FNR==NR {a[$1]; next} {if ($1=="File") {print}; if ($1 in a) {print}}' ${OUTDIR}/${Dataset}.auto.passout_QC  ${OUTDIR}/${i}.sumout_QC >${OUTDIR}/${i}.sumout_QC_onlypass
done;


###########################
### E. Deidentification ###
###########################

# For those cohorts that prefer that we do not receive the original ID, we introduce a deidentification step.

############################
# Step 1: Deidentification #
############################

# a. Make a deidentification key
declare Deidentify="${OUTDIR}/${Dataset}_deidentifykey.txt"
awk -v "name=$Dataset" 'BEGIN {OFS="\t"; i=0; print "ID_deidentified" OFS "File"} {i++; print name "_" i OFS $1}' ${OUTDIR}/${List_postQC}_adj >${Deidentify}

# b. Duplicates and relative-lists
awk -v postscript=${postscript} -v prescript=${prescript} 'BEGIN {OFS="\t"} {gsub(/.+\//,"\/",$2); gsub(/\//,"",$2); gsub(postscript,"",$2); gsub(prescript,"",$2); print}' ${Deidentify} | awk 'BEGIN {OFS="\t"} FNR==NR {a[$2]=$1; next} {if ($1 in a) {$1=a[$1]}; if ($2 in a) {$2=a[$2]}; print $0}' - ${RelativeFile}
 >${OUTDIR}/${Dataset}_DupsRelatives_key.txt

# c. Sexfile
awk -v postscript=${postscript} -v prescript=${prescript} 'BEGIN {OFS="\t"} {gsub(/.+\//,"\/",$2); gsub(/\//,"",$2); gsub(postscript,"",$2); gsub(prescript,"",$2); print}' ${Deidentify} | awk 'BEGIN {OFS="\t"} FNR==NR {a[$2]=$1; next} {if ($1 in a) {$1=a[$1]}; print $0}' - ${SexFile} >${OUTDIR}/${Dataset}_SexFile_key.txt

if [ $deidentify = "YES" ]
then
# Copy deidentified files to visualization-folder
	cp ${OUTDIR}/${Dataset}_SexFile_key.txt ${VISUALIZEDIR}/${Dataset}_SexFile.txt
	cp ${OUTDIR}/${Dataset}_DupsRelatives_key.txt ${VISUALIZEDIR}/${Dataset}_DupsRelatives.txt
else
	cp ${OUTDIR}/${Dataset}_SexFile.txt ${VISUALIZEDIR}/
	cp ${ANALYSISDIR}/DupsRelatives.txt ${VISUALIZEDIR}/${Dataset}_DupsRelatives.txt
fi
# Remaining deidentification will be done in the R-script

#############################################################
### F. Visualization of CNVs & prepare files for transfer ###
#############################################################

# Running visualization-script
sudo docker run -v ${OUTDIR}:/outdir -v ${ANALYSISDIR}:/analysisdir  -v ${VISUALIZEDIR}:/visualize -v ${ANALYSISDIR}:/pfb -v ${LRRBAFDIR}:/lrrbafdir bayramalex/enigma-cnv:latest Rscript /analysisdir/ENIGMA-CNV_visualize_v1.R ${Dataset} ${PFBname} ${CNVofInterestFile} ${Overlapref} ${OverlapMin} ${OverlapMax} ${Chip} ${deidentify}


########################
### G. Data transfer ###
########################

########################################################################
## Step 1: Putting together checklist-file including on visualization ##
########################################################################

# stop the clock
END=$(date +%s)
declare TIME=`echo $((END-START)) | awk '{print int($1/3600)":"int(($1%3600)/60)":"int($1%60)}'`

# Calculate numbers for dataset
declare IDsVisualized=`awk 'NR==1 {for (i=1; i<=NF;i++) {f[$i] = i}} { print $(f["ID_deidentified"]) }' ${ANALYSISDIR}/${Dataset}_visualize/${Dataset}_CNVcarriers_precuration_${Overlapref}.txt | sort | uniq -c | wc -l`

# a. Initiate file
rm -f ${ANALYSISDIR}/${Dataset}_visualize/${Dataset}_checklist.txt # remove file if already created
touch ${ANALYSISDIR}/${Dataset}_visualize/${Dataset}_checklist.txt
declare checklist="${ANALYSISDIR}/${Dataset}_visualize/${Dataset}_checklist.txt"

# b. Values into file

{
# Put in values for dataset etc...
echo -e "Date\t${date}
Name_of_Dataset\t${Dataset}
ResponsiblePI\t${ResponsiblePI}
Email,PI\t${Email_PI}
AnalystName\t${Analyst}
Email_Analyst\t${Email_Analyst}
TimeforAnalysis\t${TIME}
Deidentification\t${deidentify}"

#  Info for calling
awk 'END {print "Individuals_preQC\t" NR}' ${OUTDIR}/${List_preQC}
awk 'END {print "Individuals_postQC\t" NR}' ${OUTDIR}/${List_postQC}

# Generel info
echo -e "Chip_name\t$Chip
Version_if_applicable\t$Chipversion
HMM-file_used\t${HMMname}
PFB-file_used\t${PFBname}
Gcmodel-file_used\t${GCname}
No_of_individuals_with_gender_missing_in_sexfile\t$gendermissing"
echo ""

### Autosomal
awk 'END {print "Raw_dataset_Individuals_auto\t" NR-1}' ${OUTDIR}/${Dataset}.auto.sumout # individuals in the raw dataset
awk 'END {print "Filtered_dataset_Individuals_auto\t" NR}' ${OUTDIR}/${Dataset}.auto.passout # individuals in the filtered dataset
awk 'END {print "QCed_dataset_Individuals_auto\t" NR}' ${OUTDIR}/${Dataset}.auto.passout_QC # individuals in the QCed dataset
awk 'END {print "Raw_dataset_CNVs_auto\t" NR}' ${OUTDIR}/${Dataset}.auto.raw # CNVs in the raw dataset
awk 'END {print "Filtered_dataset_CNVs_auto\t" NR}' ${OUTDIR}/${Dataset}.auto.flr # CNVs in the filtered dataset
awk 'END {print "Filtered_and_merged_dataset_CNVs_auto\t" NR}' ${OUTDIR}/${Dataset}.auto.flr_mrg_final  # CNVs in the filtered AND merged dataset
awk 'END {print "Filtered_merged_and_removed_spurious_CNVs_dataset_CNVs_auto\t" NR}'  ${OUTDIR}/${Dataset}.auto.flr_mrg_spur  # CNVs in the filtered and merged and removal of spurious CNVs dataset
awk 'END {print "QCed_dataset_CNVs_auto\t" NR}' ${OUTDIR}/${Dataset}.auto.flr_QC  # CNVs in qc’ed dataset
echo ""

### X

awk 'END {print "Raw_dataset_Individuals_X\t" NR-1}' ${OUTDIR}/${Dataset}.X.sumout # Individuals in raw dataset
awk 'END {print "Filtered_dataset_Individuals_X\t" NR}' ${OUTDIR}/${Dataset}.X.passout # Individuals in the filtered dataset
awk 'END {print "Removing_individuals_based_on_autosomal_QC_Individuals_X\t" NR-1}' ${OUTDIR}/${Dataset}.X.sumout_QC_onlypass # Individuals in X after removing individuals removed in autosomal QC
awk 'END {print "Raw_dataset_CNVs_X\t" NR}' ${OUTDIR}/${Dataset}.X.raw # CNVs in the raw dataset
awk 'END {print "Filtered_dataset_CNVs_X\t" NR}' ${OUTDIR}/${Dataset}.X.flr  # CNVs in the filtered dataset
awk 'END {print "Filtered_and_merged_dataset_CNVs_X\t" NR}' ${OUTDIR}/${Dataset}.X.flr_mrg_final # CNVs in the filtered and merged dataset
awk 'END {print "Filtered_merged_and_removed_spurious_CNVs_dataset_CNVs_X\t" NR}' ${OUTDIR}/${Dataset}.X.flr_mrg_spur # CNVs in the filtered and merged and removal of spurious CNVs dataset
awk 'END {print "Removing_individuals_based_on_autosomal_QC_CNVs_X\t" NR}' ${OUTDIR}/${Dataset}.X.flr_mrg_spur_onlypass # CNVs in X after removing individuals removed in autosomal QC
awk 'END {print "QCed_dataset_CNVs_X\t" NR}' ${OUTDIR}/${Dataset}.X.flr_QC # CNVs in qc’ed dataset
echo "Please check that the above information/numbers make sense [e.g. the numbers of CNVs and/or individuals do not increase during QC]"
#### CHECKPOINTS

# Removal of nonpassed autosomal individuals
awk 'END {print "Individuals_sumoutlists_nonpassedindividualsremoved_X\t" NR-1}' ${OUTDIR}/${Dataset}.X.sumout_QC_onlypass
echo "The output should correspond to the number of individuals in the filtered autosomal dataset"

## Visualization-info

echo -e "CNVofInterestFile\t${CNVofInterestFile}
OverlapRef\t${Overlapref}
OverlapMin\t${OverlapMin}
OverlapMax\t${OverlapMax}
IDsVisualized\t${IDsVisualized}"

} >${checklist}

############################################
# Step 2: Tar the folder for file-transfer #
############################################

# All folders in the visualization folder will be transferred.

# This includes:

## Files:
   # ${Dataset}_autoX.flr_QC [Filtered CNV-dataset]
   # ${Dataset}_QC.txt [summary statistics for your raw and filtered CNV-files]
   # ${Dataset}_checklist.txt [checklist for submission]
   # ${Dataset}_CNVcarriers_precuration_0.3.txt [all CNVsofInterests that were visualized]
   # ${Dataset}_OverlapRef_0.3_TotalCNVsofInterest_precuration.txt [count of all CNVsofInterests - for overview]
   # ${Dataset}_DupsRelatives.txt [List of dups and relatives]
   # ${Dataset}_SexFile.txt [List of individuals including sex-information]
## The folders with visualization-plots:
   # singleplots_50kb
   # singleplots_cnvofinterest
   # stackplots

# Zip the visualization-folder for transfer
rm ${ANALYSISDIR}/${Dataset}_visualize.tar.gz # remove previously zipped folder
tar -zcvf ${ANALYSISDIR}/${Dataset}_visualize.tar.gz ${VISUALIZEDIR}
```

#### * Depois de 9382420984 tentativas deu certo! Comentei a parte do PFB e GCmodel porque já haviam sido gerados em uma das rodadas. 

## **Outros arquivos para mandar**

### **D. Covariate files**
The covariate files are produced by manually constructing a covariate-file and obtaining ancestry covariates following the standard protocols of ENIGMA genetics. 
Please create a csv-file called ${Dataset]_Covar_ENIGMACNV.csv
Columns: SubjID	GeneticID	PID	     MID	DiseaseType	  AffectionStatus	Affectionstatus2	Age	  Sex	ScannerSite
		  Subj1	  A2035	   A1003	A6008	Population			0			Healthy_bipolar		30.0  male     Oslo2

General NOTES on the columns
•     Missing values should be coded as NA. 
•     Please include individuals, even if they are missing one or two covariates (say AffectionStatus) – the data might still be valuable. Just fill out the missing values as NA.

Specific NOTES on the columns
•     SubjID: Must match the SubjID in the imaging files.
• GeneticID: The ID used in the ENIGMA-CNV-calling protocol (“ENIGMA- CNV_CNVCalling_Protocol_final.sh”) prior to de-identification. If this is the same as the SubjID, fill in the column accordingly.
•     Sex: Must be coded as follows: Male=”male”, Female= ”female”. Note sex-info is also needed for the CNV-calling above.
•     DiseaseType: Type of cohort. 
ENIGMA-CNV receives both population cohorts (or volunteer-cohorts) and disease-cohorts (e.g. epilepsy, psychiatric disease, dementia).
•        For population-studies, please note the best fitting term (copy to all fields):
o Twinstudy
o Population
o Volunteers
o Familystudy
•        For case-controls-studies, please note the best fitting term (copy to all fields): o E.g bipolar, schizophrenia, ADHD, autism, dementia, stroke, psychosis, epilepsy… ## Colocaremos ANX, MDD, BD, PTSD, ADHD?
•     AffectionStatus:  a binary indicator where Patient = 1 & Control = 0. If your cohort does not have patients, code = 0. 
•        AffectionStatus2
•        For population-studies: please code as NA in all columns
•        For case-controls-studies:
o Please code patients with
§        a standardized version of the diagnosis per individual such as DSM or ICD (We are particularly interested in ‘Mental, behavioural and neurodevelopmental disorders’, corresponding to ICD-10 F01-99 and ‘Diseases of the nervous system, corresponding to ICD-10 G00-G99.)
§        If lacking such detailed information, SCZ (schizophrenia), BD (bipolar disorder), ADHD, dementia, stroke or the likes will suffice.
o Please code controls as:
§ Healthy_psych, Healthy_epilepsy, Healthy_dementia (controls have been ascertained/screened for lack of e.g. psychiatric diseases, epilepsy or dementia)
•     ScannerSite: The name of your scanner site.
OBS: Make sure the ScannerSite name in Covar_ENIGMACNV.csv correlate with that/those put into n ScannerInfoSheet.xlxs (see sMRI protocol).

### **MDS-covariates**
MDS-covariate-file: HM3_b37mds2R.mds.csv (version July 13, 2017), alternatively HM3mds2R.mds.csv (version July 27, 2012).
This file contains covariates for the ancestry of each subject in your cohort.
If you previously participated in ENIGMA-genetics (and followed the protocols of ENIGMA), you have this file already.
If not, use this link:
https://enigma.ini.usc.edu/wp-content/uploads/2020/02/ENIGMA-1KGP_p3v5- Cookbook_20170713.pdf

OBS: To obtain MDS-covariates, you only need to run the “Multi-dimensional Scaling (MDS) Protocol” part of the code in the ENIGMA imputation protocol, NOT the imputation part. 

The output file, HM3_b37mds2R.mds.csv, is a spreadsheet containing the following columns:
Family ID (FID), individual ID (IID), 4 MDS components (C1, C2, C3 and C4), and PLINK’s assigned solution code (SOL). If you have more MDS components available, please provide us with up to 20.

```
wget "http://genepi.qimr.edu.au/staff/sarahMe/enigma/MDS/HM3_b37.bed.gz" 
wget "http://genepi.qimr.edu.au/staff/sarahMe/enigma/MDS/HM3_b37.bim.gz" 
wget "http://genepi.qimr.edu.au/staff/sarahMe/enigma/MDS/HM3_b37.fam.gz" 
```
```
plink --bfile BHRC_Probands --hwe 1e-6 --geno 0.05 --maf 0.01 --noweb --make-bed --out BHRC_Probands_filtered
```












### **E. Return data to ENIGMA-CNV**
Files to transfer:

BHRCS_visualize.tar.gz [CNV-calling data]

BHRCS_sMRI.tar.gz [sMRI calling data]

HM3_b37mds2R.mds.csv [MDS-data]

BHRCS_Covar_ENIGMACNV.csv [covariates]

