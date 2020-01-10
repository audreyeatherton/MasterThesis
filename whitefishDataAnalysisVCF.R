# Work on the vcf
# VCF too big for Lappie so done on axiom
# Also Lappie doesn't have hierfstat


#######################
# Plot Heterozygosity #
#######################

# Run from cluster
library(hierfstat)
library(vcfR)
library(Matrix)
library(lme4)
library(ggplot2)
library(plyr)
library(dplyr)
library(tidyr)
library(gplots)
library(calibrate)

# Set things up
rm(list=ls())
getwd()
setwd("/home/emma/Desktop/Project/FinalData/")

# Load the genetic data

whitefish <- read.vcfR('whitefish.vcf')

# refgenome <- ape::read.dna('whitefish_genome.fasta', format = 'fasta')
# 
# chrom <- create.chromR(name = "chrom", vcf = whitefish, seq = refgenome)


############################
# Plot SNPs per Chromosome #
############################






# Extract genotypes and ID of the fish
gen <- extract.gt(whitefish, 'GT')
gen <- extract.gt(chrom, 'GT')

# Make it a dataframe
gen1 <- as.data.frame(gen)

# Transpose to make compatible with hierfstat
gen2 <- t(gen1)
# write.csv(gen2, file = "sireGenotypes.csv")
# gen2 <- read.csv('sireGenotypes.csv')
# Modify the genotype value to homozygous = 0 or heterozygous = 1
gen1113 <- mapvalues(gen2$X1113, from = c("0/0", "0/1", "1/1"), to = c("0", "1", "0"))

gen3 <- as.data.frame(cbind(gen1102, gen1103, gen1105, gen1106, gen1107, gen1108, gen1109,
                            gen1110, gen1111, gen1112, gen1113, gen1114, gen1115, gen1117))
# write.csv(gen3, file = "sireGenotypes0s1s.csv")
gen3 <- read.csv("sireGenotypes0s1s.csv")

# Calculate the allele frequency for each SNP
maf <- colMeans(gen3, na.rm = T)/2
# Make the max maf 0.5, so half of the HW equilibrium
maf[maf>0.5] <- 1 - maf[maf>0.5]

# Estimate the heterozygosity per sample
ho <- apply(gen3, 2, function(x)sum(x==1, na.rm = T)/length(which(!is.na(x))))
# Plot observed heterozygosity and allele frequency, and the expected values in HW equilibrium
plot(maf,ho)
x <- seq(0, 0.5, 0.001)
lines(x, 2*x*(1-x))
hist(maf)

hetdta <- as.data.frame(cbind(ho, maf))
ggplot(data = hetdta, aes(maf, ho)) +
  geom_point() 
