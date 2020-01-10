# Pipeline used for preprocessing of whitefish whole genome sequencing data in March to September 2019
# Relies on already publicly available reference genome (which came out in August 2019)





# Part 1   Get raw data
# You can't analyze data you don't have





# Download data

# Rename files


# Move raw data from old cluster to new cluster
# movedat -i aatherto@dexvital.unil.ch:/vital-it/archive/dee/wedekind/aatherto/aatherto/fasta/ /scratch/wally/FAC/FBM/DEE/cwedekin/whitefishhallwil/aatherto/sequencing



# Slurm header
#!/bin/sh

#SBATCH --account cwedekin_whitefishhallwil
#SBATCH --mail-user audrey.atherton@unil.ch
#SBATCH --mail-type ALL
#SBATCH --partition ax-normal
#SBATCH --time 24:00:00
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 8
#SBATCH --mem 16G
#SBATCH --job-name filter_whitefish
#SBATCH --export NONE

# Load modules
module load Bioinformatics/Software/vital-it
module add UHTS/Quality_control/fastqc/0.11.7;
module add UHTS/Aligner/bwa/0.7.17;
module add UHTS/Analysis/samtools/1.8;
module add UHTS/Analysis/picard-tools/2.18.11;
module add UHTS/Analysis/GenomeAnalysisTK/4.1.3.0;
module add UHTS/Analysis/vcftools/0.1.15;
module add UHTS/Analysis/plink/1.90;
module add UHTS/Analysis/stacks/1.48;
module add UHTS/Analysis/EPACTS/3.2.6;
module add R;




# Part 2   Data preprocessing
# Make sure the data you analyze is of high enough quality to mean more than nothing



# Quality control on the raw data
ls -v1 | grep .fastq.gz | xargs -d'\n' -n 1 fastqc -t 2
# list all files in the folder
# take only the .fastq.gz
# take only one at a time as an argument for fastqc
# -o set output directory
# fastqc option -t is for the number of files that can be processed simultaneously (2)



# Align raw data and turn to bam
bwa mem -t 20 /scratch/wally/FAC/FBM/DEE/cwedekin/whitefishhallwil/aatherto/sequencing/whitefish_genome.fasta /scratch/wally/FAC/FBM/DEE/cwedekin/whitefishhallwil/aatherto/sequencing/yoda10/1102_roh_L1_R1.fastq /scratch/wally/FAC/FBM/DEE/cwedekin/whitefishhallwil/aatherto/sequencing/yoda10/1102_roh_L1_R2.fastq | samtools view -Sb - > /scratch/wally/FAC/FBM/DEE/cwedekin/whitefishhallwil/aatherto/sequencing/yoda101102_L1.bam
# Format
# bwa mem -t (number of cores) /path/to/reference/genome /path/to/first/run/of/sequence/file.fastq /path/to/second/run/of/sequence/file.fastq | samtools view -Sb - > /path/to/output/file.bam
# Piped to avoid the intermediary .sam file, which is bigger and which we don't work on directly



# Clean with picard tools
picard-tools CleanSam I=/scratch/wally/FAC/FBM/DEE/cwedekin/whitefishhallwil/aatherto/sequencing/yoda10/bam/1102_L1.bam O=/scratch/wally/FAC/FBM/DEE/cwedekin/whitefishhallwil/aatherto/sequencing/yoda10/cleanedBam/1102_L1_cleaned.bam
# Format
# picard-tools CleanSam I=/path/to/input/file.bam O=/path/to/output/file.bam



# Add read groups
# Python function to set the read groups 
def picardAddReadGroups(runname):
  # Load libraries
  import os

  # Make loop
  for file in os.listdir():
    # Import data from bam files
    path = os.path.abspath(file)
    filename = file

    # Build necessary things
    lanenumber = file[file.find("L") + 1]
      # Will have 1 for yoda10 and 9 for yoda9
      # If working with double digits, change to following
    # lanenumber = file[file.find("L") + 1:file.find("L") + 2]
    runnumber = path[path.find(runname) + len(runname)]

    # Name variables for picard
    i = path
	# Path to input file
    o = "../../cleanRGBam/" + file
	# Path to output file
    rgid = lanenumber + runnumber
	# Name of run group (here lane and yoda number)
    rglb = "lib1"
	# Read group library default, not sure what it means
    rgpl = "illumina"
	# Platform used for sequencing
    rgpu = runname
	# Name of run (here yoda)
    rgsm = file[0:4]
	# Name of individual subject (here between 1102 and 1117)

    # Run picard
    pic = ("picard-tools AddOrReplaceReadGroups I=%s O=%s RGID=%s RGLB=%s RGPL=%s RGPU=%s RGSM=%s" % (i, o, rgid, rglb, rgpl, rgpu, rgsm))
    os.system(pic)

picardAddReadGroups("yoda")
# Call function from within the same .py file

# Call of the python code
python3 addReadGroup.py



# Fix paired ends
samtools fixmate /scratch/wally/FAC/FBM/DEE/cwedekin/whitefishhallwil/aatherto/sequencing/yoda10/cleanRGBam/1102_L1_cleaned.bam /scratch/wally/FAC/FBM/DEE/cwedekin/whitefishhallwil/aatherto/sequencing/yoda10/fixedPEBam/1102_L1_cleaned_RG_fixPE.bam
# Format
# samtools fixmate /path/to/input/file.bam /path/to/output/file.bam



# Remove secondary alignments
samtools view -bh -F 256 /scratch/wally/FAC/FBM/DEE/cwedekin/whitefishhallwil/aatherto/sequencing/yoda10/fixedPEBam/1102_L1_cleaned_RG_fixPE.bam -o /scratch/wally/FAC/FBM/DEE/cwedekin/whitefishhallwil/aatherto/sequencing/yoda10/secRemBam/1102_L1_cleaned_RG_fixedPE_noSec.bam
# Format
# samtools view -bh -F 256 /path/to/input/file.bam -o /path/to/output/file.bam



# Sort resulting bam files
samtools sort /scratch/wally/FAC/FBM/DEE/cwedekin/whitefishhallwil/aatherto/sequencing/yoda10/secRemBam/1102_L1_cleaned_RG_fixedPE_noSec.bam -o /scratch/wally/FAC/FBM/DEE/cwedekin/whitefishhallwil/aatherto/sequencing/yoda10/sortBam/1102_L1_cleaned_RG_fixedPE_noSec_sorted.bam
# Format
# samtools sort /path/to/input/file.bam -o /path/to/output/file.bam
# Needs reference genome in the same directory to run



# Mark duplicates
picard-tools MarkDuplicates I=/scratch/wally/FAC/FBM/DEE/cwedekin/whitefishhallwil/aatherto/sequencing/yoda10/sortBam/1102_L1_cleaned_RG_fixedPE_noSec_sorted.bam O=/scratch/wally/FAC/FBM/DEE/cwedekin/whitefishhallwil/aatherto/sequencing/yoda10/markedDupBam/1102_L1_cleaned_RG_fixedPE_noSec_sorted_markedDup.bam M=/scratch/wally/FAC/FBM/DEE/cwedekin/whitefishhallwil/aatherto/sequencing/yoda10/markedDupBam/1102_L1_cleaned_RG_fixedPE_noSec_sorted_markedDup_metrics.txt
# Format
# picard-tools MarkDuplicates I=/path/to/input/file.bam O=/path/to/output/file.bam M=/path/to/metrics/output/file.txt



# GATK referencing of genome:
picard-tools CreateSequenceDictionary R=/scratch/wally/FAC/FBM/DEE/cwedekin/whitefishhallwil/aatherto/sequencing/whitefish_genome.fasta O=/scratch/wally/FAC/FBM/DEE/cwedekin/whitefishhallwil/aatherto/sequencing/whitefish_genome.dict
# Format
# picard-tools CreateSequenceDictionary R=/path/to/reference/genome/input.fasta O=/path/to/reference/genome/output.dict
# Make sure all the genome reference files are in the same directory

# Make indexed fasta file (necessary for reordering with picard-tools!)
samtools faidx /scratch/wally/FAC/FBM/DEE/cwedekin/whitefishhallwil/aatherto/sequencing/whitefish_genome.fasta
# Format
# samtools faidx /path/to/reference/genome/input.fasta
# Output is an automatically generated .fai file necessary for reordering of bam file with picard-tools


# Reorder bam along referenced genome (GATK)
picard-tools ReorderSam I=/scratch/wally/FAC/FBM/DEE/cwedekin/whitefishhallwil/aatherto/sequencing/yoda10/markedDupBam/1102_L1_cleaned_RG_fixedPE_noSec_sorted_markedDup.bam O=/scratch/wally/FAC/FBM/DEE/cwedekin/whitefishhallwil/aatherto/sequencing/yoda10/reorderedBam/1102_L1_cleaned_RG_fixedPE_noSec_sorted_markedDup_reordered.bam R=whitefish_genome.fasta
# Format
# picard-tools ReorderSam I=/path/to/input/file.bam O=/path/to/output/file.bam R=/path/to/reference/genome.fasta



# Generate indexed bam files, bai
samtools index /scratch/wally/FAC/FBM/DEE/cwedekin/whitefishhallwil/aatherto/sequencing/yoda10/reorderedBam/1102_L1_cleaned_RG_fixedPE_noSec_sorted_markedDup_reordered.bam
# Format
# samtools index /path/to/input/file.bam
# Needs reference genome in same directory, output is a .bai file



# Validate bam
picard-tools ValidateSamFile I=/scratch/wally/FAC/FBM/DEE/cwedekin/whitefishhallwil/aatherto/sequencing/yoda10/reorderedBam/1102_L1_cleaned_RG_fixedPE_noSec_sorted_markedDup_reordered.bam O=/scratch/wally/FAC/FBM/DEE/cwedekin/whitefishhallwil/aatherto/sequencing/yoda10/validBam/1102_L1_cleaned_RG_fixedPE_noSec_sorted_markedDup_reordered_valid.txt IGNORE=INVALID_FLAG_MATE_UNMAPPED IGNORE= MISMATCH_FLAG_MATE_UNMAPPED IGNORE=MISMATCH_MATE_CIGAR_STRING IGNORE= MATE_CIGAR_STRING_INVALID_PRESENCE
# Format
# picard-tools ValidateSamFile I=/path/to/input/file.bam O=/path/to/output/file.txt IGNORE=INVALID_FLAG_MATE_UNMAPPED IGNORE= MISMATCH_FLAG_MATE_UNMAPPED IGNORE=MISMATCH_MATE_CIGAR_STRING IGNORE= MATE_CIGAR_STRING_INVALID_PRESENCE
# Had only one error about the mate cigar string being unmapped
# Ignored the error because it was expected and we wanted to know if any other errors existed



# Merge lanes bam
picard-tools MergeSamFiles I=/scratch/axiom/FAC/FBM/DEE/cwedekin/whitefishhallwil/aatherto/bamEndPrep/yoda10/1102_L1_cleaned_RG_fixedPE_noSec_sorted_markedDup_reordered.bam I=/scratch/axiom/FAC/FBM/DEE/cwedekin/whitefishhallwil/aatherto/bamEndPrep/yoda10/1102_L2_cleaned_RG_fixedPE_noSec_sorted_markedDup_reordered.bam I=/scratch/axiom/FAC/FBM/DEE/cwedekin/whitefishhallwil/aatherto/bamEndPrep/yoda10/1102_L3_cleaned_RG_fixedPE_noSec_sorted_markedDup_reordered.bam I=/scratch/axiom/FAC/FBM/DEE/cwedekin/whitefishhallwil/aatherto/bamEndPrep/yoda10/1102_L4_cleaned_RG_fixedPE_noSec_sorted_markedDup_reordered.bam I=/scratch/axiom/FAC/FBM/DEE/cwedekin/whitefishhallwil/aatherto/bamEndPrep/yoda10/1102_L5_cleaned_RG_fixedPE_noSec_sorted_markedDup_reordered.bam I=/scratch/axiom/FAC/FBM/DEE/cwedekin/whitefishhallwil/aatherto/bamEndPrep/yoda10/1102_L6_cleaned_RG_fixedPE_noSec_sorted_markedDup_reordered.bam I=/scratch/axiom/FAC/FBM/DEE/cwedekin/whitefishhallwil/aatherto/bamEndPrep/yoda10/1102_L7_cleaned_RG_fixedPE_noSec_sorted_markedDup_reordered.bam I=/scratch/axiom/FAC/FBM/DEE/cwedekin/whitefishhallwil/aatherto/bamEndPrep/yoda10/1102_L8_cleaned_RG_fixedPE_noSec_sorted_markedDup_reordered.bam I=/scratch/axiom/FAC/FBM/DEE/cwedekin/whitefishhallwil/aatherto/bamEndPrep/yoda9/1102_L1_cleaned_RG_fixedPE_noSec_sorted_markedDup_reordered.bam I=/scratch/axiom/FAC/FBM/DEE/cwedekin/whitefishhallwil/aatherto/bamEndPrep/yoda9/1102_L2_cleaned_RG_fixedPE_noSec_sorted_markedDup_reordered.bam I=/scratch/axiom/FAC/FBM/DEE/cwedekin/whitefishhallwil/aatherto/bamEndPrep/yoda9/1102_L3_cleaned_RG_fixedPE_noSec_sorted_markedDup_reordered.bam I=/scratch/axiom/FAC/FBM/DEE/cwedekin/whitefishhallwil/aatherto/bamEndPrep/yoda9/1102_L4_cleaned_RG_fixedPE_noSec_sorted_markedDup_reordered.bam I=/scratch/axiom/FAC/FBM/DEE/cwedekin/whitefishhallwil/aatherto/bamEndPrep/yoda9/1102_L5_cleaned_RG_fixedPE_noSec_sorted_markedDup_reordered.bam I=/scratch/axiom/FAC/FBM/DEE/cwedekin/whitefishhallwil/aatherto/bamEndPrep/yoda9/1102_L6_cleaned_RG_fixedPE_noSec_sorted_markedDup_reordered.bam I=/scratch/axiom/FAC/FBM/DEE/cwedekin/whitefishhallwil/aatherto/bamEndPrep/yoda9/1102_L7_cleaned_RG_fixedPE_noSec_sorted_markedDup_reordered.bam I=/scratch/axiom/FAC/FBM/DEE/cwedekin/whitefishhallwil/aatherto/bamEndPrep/yoda9/1102_L8_cleaned_RG_fixedPE_noSec_sorted_markedDup_reordered.bam O=/scratch/axiom/FAC/FBM/DEE/cwedekin/whitefishhallwil/aatherto/bamEndPrep/1102_end_preprocessing.bam CREATE_INDEX=true
# picard-tools MergeSamFiles I=/path/to/ind1/lane1.bam I=I=/path/to/ind1/lane2.bam [...] I=/path/to/ind1/lanen.bam O=/path/to/output/file.bam CREATE_INDEX=true



# Part 3   SNP Calling
# Focus on the parts that are different from the reference genome
# Smaller files with only the relevant information (I hope)
	# Note : One can also merge bam files together before SNP calling
	# Using picard-tools MergeSamFiles
	# Code for this in file mergeBamPreSNPCalling.sh




# Make gVCF files with haplotype caller (1 per lane, per individual)
GenomeAnalysisTK HaplotypeCaller -R /scratch/wally/FAC/FBM/DEE/cwedekin/whitefishhallwil/aatherto/sequencing/whitefish_genome.fasta -I /scratch/wally/FAC/FBM/DEE/cwedekin/whitefishhallwil/aatherto/sequencing/yoda10/reorderedBam/1102_L1_cleaned_RG_fixedPE_noSec_sorted_markedDup_reordered.bam -L /scratch/axiom/FAC/FBM/DEE/cwedekin/whitefishhallwil/aatherto/bamEndPrep/chrom01to04.intervals --emit-ref-confidence GVCF --genotyping-mode DISCOVERY -O /scratch/wally/FAC/FBM/DEE/cwedekin/whitefishhallwil/aatherto/sequencing/yoda10/gvcf/1102_L1.g.vcf
# Format
# GenomeAnalysisTK HaplotypeCaller -R /path/to/reference/genome.fasta -I /path/to/input/file.bam -L /path/to/chromosome/interval/file.intervals --emit-ref-confidence GVCF --genotyping-mode DISCOVERY -O /path/to/output/file.g.vcf
# --emit-ref-confidence GVCF and --genotyping-mode DISCOVERY are VITAL! Without them, GATK won't recognize the output!




# Merge GVCF files 
GenomeAnalysisTK CombineGVCFs -R /scratch/axiom/FAC/FBM/DEE/cwedekin/whitefishhallwil/aatherto/genome/whitefish_genome.fasta -V /scratch/axiom/FAC/FBM/DEE/cwedekin/whitefishhallwil/aatherto/gvcfs/gvcfMergedBamSepChrom/1102_chrom01to04_new.g.vcf -V /scratch/axiom/FAC/FBM/DEE/cwedekin/whitefishhallwil/aatherto/gvcfs/gvcfMergedBamSepChrom/1102_chrom05to08_new.g.vcf -V /scratch/axiom/FAC/FBM/DEE/cwedekin/whitefishhallwil/aatherto/gvcfs/gvcfMergedBamSepChrom/1102_chrom09to12_new.g.vcf -V /scratch/axiom/FAC/FBM/DEE/cwedekin/whitefishhallwil/aatherto/gvcfs/gvcfMergedBamSepChrom/1102_chrom13to16_new.g.vcf -V /scratch/axiom/FAC/FBM/DEE/cwedekin/whitefishhallwil/aatherto/gvcfs/gvcfMergedBamSepChrom/1102_chrom17to20_new.g.vcf -V /scratch/axiom/FAC/FBM/DEE/cwedekin/whitefishhallwil/aatherto/gvcfs/gvcfMergedBamSepChrom/1102_chrom21to24_new.g.vcf -V /scratch/axiom/FAC/FBM/DEE/cwedekin/whitefishhallwil/aatherto/gvcfs/gvcfMergedBamSepChrom/1102_chrom25to28_new.g.vcf -V /scratch/axiom/FAC/FBM/DEE/cwedekin/whitefishhallwil/aatherto/gvcfs/gvcfMergedBamSepChrom/1102_chrom29to32_new.g.vcf -V /scratch/axiom/FAC/FBM/DEE/cwedekin/whitefishhallwil/aatherto/gvcfs/gvcfMergedBamSepChrom/1102_chrom33to36_new.g.vcf -V /scratch/axiom/FAC/FBM/DEE/cwedekin/whitefishhallwil/aatherto/gvcfs/gvcfMergedBamSepChrom/1102_chrom37to40_new.g.vcf -O /scratch/axiom/FAC/FBM/DEE/cwedekin/whitefishhallwil/aatherto/gvcfs/gvcfMergedBamSepToInd/1102_long.g.vcf
# Way it was done
	# Merged all chromosomes of a single individual
	# Merged the first 6 individuals, and separately the last 8 individuals (couldn't run all 14 at once)
	# Merge the 2 final files into a single GVCF file containing every lane of every individual for both yoda9 and yoda10
# Format
# GenomeAnalysisTK CombineGVCFs -R /path/to/reference/genome.fasta -V /path/to/input/file1.g.vcf -V /path/to/input/file2.g.vcf -V /path/to/input/file3.g.vcf -O /path/to/output/file.g.vcf



# Turn merged GVCF to VCF
GenomeAnalysisTK GenotypeGVCFs -R /scratch/axiom/FAC/FBM/DEE/cwedekin/whitefishhallwil/aatherto/whitefish_genome.fasta -V /scratch/axiom/FAC/FBM/DEE/cwedekin/whitefishhallwil/aatherto/whitefishLongQ.g.vcf -O /scratch/axiom/FAC/FBM/DEE/cwedekin/whitefishhallwil/aatherto/whitefishLong.vcf
# Format
# GenomeAnalysisTK GenotypeGVCFs -R /path/to/reference/genome.fasta -V /path/to/input/file.g.vcf -O /path/to/output/file.vcf
# Only takes one input file





# Part 4   Filtering
# Keep only what is likely to be true variation, eliminate even more potential errors





# Hard Filtering with GATK's VariantFiltration
GenomeAnalysisTK VariantFiltration -V /scratch/axiom/FAC/FBM/DEE/cwedekin/whitefishhallwil/aatherto/vcfAllInOne/whitefishLong.vcf -filter "QD<2.0" --filter-name "QD2" -filter "QUAL<30.0" --filter-name "QUAL30" -filter "SOR>3.0" --filter-name "SOR3" -filter "FS>60.0" --filter-name "FS60" -filter "MQ<40.0" --filter-name "MQ40" -filter "MQRankSum<-12.5" --filter-name "MQRankSum-12.5" -filter "ReadPosRankSum<-8.0" --filter-name "ReadPosRankSum-8" -O /scratch/axiom/FAC/FBM/DEE/cwedekin/whitefishhallwil/aatherto/filter/whitefishQualFilterRedo.vcf
# Format
# GenomeAnalysisTK VariantFiltration -V /path/to/input/file.vcf -filter "QD<2.0" --filter-name "QD2" -filter "QUAL<30.0" --filter-name "QUAL30" -filter "SOR>3.0" --filter-name "SOR3" -filter "FS>60.0" --filter-name "FS60" -filter "MQ<40.0" --filter-name "MQ40" -filter "MQRankSum<-12.5" --filter-name "MQRankSum-12.5" -filter "ReadPosRankSum<-8.0" --filter-name "ReadPosRankSum-8" -O /path/to/output/file.vcf
# Default values and filters found : https://software.broadinstitute.org/gatk/documentation/article?id=23216 and https://gatkforums.broadinstitute.org/gatk/discussion/2806/howto-apply-hard-filters-to-a-call-set
# QD		   Normalize the variant quality to avoid inflation caused by deep coverage (if coverage is deeper, the variant quality seems better than it truly is)
# QUAL		   Quality
# MQ		   Root mean square mapping quality over all the reads at the site (if the read doesn't map well, 
# FS		   Strand bias. Does forward or reverse strand have more alternate allelles (should be about equal, since they're supposedly reverse copies of each other). Penalized ends of exons
# SOR		   Strand bias. Similar to FS but doesn't penalize ends of exons, which tend to be covered only in one direction
# MQRankSum	   Compares the mapping qualities of reads supporting the reference allele and the alternate allele (negative -> reference has better mappi:wqng than alternate, so alternate may be error)
# ReadPosRankSum   Checks whether the positions of the reference and alternate alleles are different within the reads (negative -> alternate allele is at the end of reads more often, so likely an error)



# Make new vcf file excluding the variants that didn't pass the technical filter
GenomeAnalysisTK SelectVariants -R /scratch/axiom/FAC/FBM/DEE/cwedekin/whitefishhallwil/aatherto/genome/whitefish_genome.fasta -V /scratch/axiom/FAC/FBM/DEE/cwedekin/whitefishhallwil/aatherto/filter/whitefishTechFilt.vcf -O /scratch/axiom/FAC/FBM/DEE/cwedekin/whitefishhallwil/aatherto/filter/whitefishEndTechFilt.vcf --exclude-filtered
# Format
# GenomeAnalysisTK SelectVariants -R /path/to/reference/file.fasta -V /path/to/input/file.vcf -O /path/to/output/file.vcf --exclude-filtered
# Careful, the website says to use --excludeFiltered, but that won't run



# IMPORTANT! 
# Order of one-at-a-time filtering
	# --hwe 0.05		vcftools
	# --maf 0.05		vcftools
	# --min-meanDP 10	vcftools
	# --max-meanDP 60	vcftools
	# --minGQ 20		vcftools
	# --max-missing 0.9	vcftools
	# --remove-indels	vcftools
	# --max-alleles 2	vcftools
	# excess heterozygotes	stacks


# Filter  depth and missing data with VCFTools
# Filter only ONE thing at a time for this step! Boundaries will depend on the data, species, biological question, ...
vcftools --gzvcf /scratch/axiom/FAC/FBM/DEE/cwedekin/whitefishhallwil/aatherto/filterVarFiltr/whitefishFilteredQual.vcf.gz --out /scratch/axiom/FAC/FBM/DEE/cwedekin/whitefishhallwil/aatherto/filterVCFTools/whitefishFilteredQualAndVCFT --max-missing 0.3 --recode
# Format 
# vcftools --gzvcf /path/to/input/file.vcf.gz --out /path/to/output/file --option value --recode
# IMPORTANT! The --recode option is super important when filtering, without it no output file is written!



# Rename chromosomes into something understandable
sed -i "s/LR664344.1/01/g" whitefish.vcf
# All chromosomes were thus renamed according to their official number on the official genome (https://www.ncbi.nlm.nih.gov/genome/?term=Coregonus as of 2019 December 6th)



# Rename Variant IDs with bcftools
bcftools annotate --set-id +'%CHROM\_%POS'  whitefishFinal.vcf --output whitefishFinalNamed.vcf



# Filter out excess heterosity (80% max)
populations -V /scratch/axiom/FAC/FBM/DEE/cwedekin/whitefishhallwil/aatherto/filter/whitefishFinalNamed.vcf -O /scratch/axiom/FAC/FBM/DEE/cwedekin/whitefishhallwil/aatherto/stacksHet -M /scratch/axiom/FAC/FBM/DEE/cwedekin/whitefishhallwil/aatherto/filter/pop_map.txt -p 1 --max_obs_het 0.80 --vcf
# Format
# populations -V /path/to/input/file.vcf -O /path/to/output/directory -M /path/to/population/map/file.txt -p 1 --max_obs_het 0.80 --vcf
# population map is a text file containing one line per individual following the pattern "ind_ID '\t' pop_nb"
# Man page : http://catchenlab.life.illinois.edu/stacks/comp/populations.php
# Used Stacks 1.48, extremely RAM greedy! Needed 500G for it to run on a 4.6G vcf (runtime 1h15)



# Get Multi Locus Heterozygosity
vcftools --vcf whitefishFinalNamed.vcf --het --out whitefishHet



# Turn VCF to BED
plink --vcf /scratch/axiom/FAC/FBM/DEE/cwedekin/whitefishhallwil/aatherto/bed/whitefishFiltered.vcf --make-bed --chr-set 40 --out /scratch/axiom/FAC/FBM/DEE/cwedekin/whitefishhallwil/aatherto/bed/whitefishName
# Whole VCF





# Part 5	ROH calling with Plink




# Calling ROH with Plink
plink --bfile whitefishAll --homozyg --homozyg-snp 100 --homozyg-kb 100 --homozyg-density 5000 --homozyg-gap 5000 --homozyg-window-snp 100 --homozyg-window-het 10 --homozyg-window-missing 25 --homozyg-window-threshold 0.05 --chr-set 40 --out whitefishAllROH
# Format
# plink --bfile inputFileName --homozyg --homozyg-snp 100 --homozyg-kb 100 --homozyg-density 5000 --homozyg-gap 5000 --homozyg-window-snp 100 --homozyg-window-het 10 --homozyg-window-missing 25 --homozyg-window-threshold 0.05 --chr-set 40 --out outputFileName
# --homozyg 				function for ROH calling
# --homozyg-snp 100 			minimum SNPs per ROH
# --homozyg-kb 100			minimum size of ROHs, in kb
# --homozyg-density 5000		minimum SNP density (which we excluded)
# --homozyg-gap 5000			minimum inter-snp gap (which we excluded)
# --homozyg-window-snp 100		number of SNPs in scanning window 
# --homozyg-window-het 10		maximum number of heterozygous SNPs in ROH
# --homozyg-window-missing 25		maximum number of missing genotypes in ROH
# --homozyg-window-threshold 0.05	threshold for SNP to be included in ROH
# --chr-set 40				number of chromosomes in our genome (not human)
# Man page : https://www.cog-genomics.org/plink/1.9/ibd
# 21'000 ROHs, includes small to large ROHs



plink --bfile whitefishAll --homozyg --homozyg-snp 100 --homozyg-kb 100 --homozyg-density 5000 --homozyg-gap 5000 --homozyg-window-snp 100 --homozyg-window-het 2 --homozyg-window-missing 25 --homozyg-window-threshold 0.05 --chr-set 40 --out whitefishAllROHTest3

plink --bfile whitefishAllExHet05 --homozyg --homozyg-snp 100 --homozyg-kb 100 --homozyg-density 5000 --homozyg-gap 5000 --homozyg-window-snp 100 --homozyg-window-het 10 --homozyg-window-missing 25 --homozyg-window-threshold 0.05 --chr-set 40 --out whitefishAllROHExHet05

plink --bfile whitefishAllExHet05 --homozyg --homozyg-snp 100 --homozyg-kb 100 --homozyg-density 5000 --homozyg-gap 5000 --homozyg-window-snp 100 --homozyg-window-het 2 --homozyg-window-missing 25 --homozyg-window-threshold 0.05 --chr-set 40 --out whitefishAllROHExHet05Strict
