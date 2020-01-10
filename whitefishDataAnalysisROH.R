# Whitefish data analysis with Runs of Homozygosity
# Whole genome data
# Audrey Atherton (Based on code from Christian DeGuttry)
# 2018-2019


##############
### Set-Up ###
##############

# Load necessary libraries
library(vcfR)
library(Matrix)
library(lme4)
library(ggplot2)
library(viridis)
library(hrbrthemes)
library(plyr)
library(dplyr)
library(tidyr)
library(calibrate)
library(tidyverse)
library(RColorBrewer)
library(outliers)


# Clear environment
rm(list=ls())
# Set working directory
setwd("/home/emma/Desktop/Project/FinalData/")


# Load data 
# Offspring data
dta <- read.csv('Pf2_20_12_all measurements_23_10_19.csv', header = T)
# Exclude the individuals from Sires from whom we don't have whole genome data
dta <- subset(dta, Sire!=c("1101"))
dta <- subset(dta, Sire!=c("1104"))
dta <- subset(dta, Sire!=c("1116"))
dta <- subset(dta, Sire!=c("1118"))

# Sire data
dads <- read.csv('20_12 Breeders final (09_08) (1).csv', header = T)
# Exclude the individuals from Sires from whom we don't have whole genome data
dads <- subset(dads, Sire!=(""))
dads <- subset(dads, Sire!=c("1101"))
dads <- subset(dads, Sire!=c("1104"))
dads <- subset(dads, Sire!=c("1116"))
dads <- subset(dads, Sire!=c("1118"))

# ROH data
ROH <- read.table('whitefishAllROHExHet05.txt', sep = '\t', header = T) 




#############################
# Graph SNPs per Chromosome #
#############################

Chr <- seq(1, 40, 1)
ChrSizeMb <- c(93.46, 43.33, 42.76, 92.22, 45.59,   43.69, 65.39, 68.14, 65.19, 60.47,
               63.18, 63.88, 57.74, 47.26, 55.64,   54.03, 54.22, 51.95, 59.91, 54.34,
               52.95, 56.86, 52.02, 50.33, 51.03,   50.92, 48.68, 46.67, 48.98, 48.68,
               48.45, 44.62, 44.66, 40.73, 42.61,   42.01, 43.66, 36.77, 33.96, 01.09)
# Data from ncbi genome  https://www.ncbi.nlm.nih.gov/genome/?term=Coregonus
SNPperChr <- c(385294, 207515, 198268, 417540, 217326,   187133, 248269, 329502, 269669, 307815,
               275454, 282970, 258287, 243851, 236867,   241301, 213398, 237778, 277330, 196574,
               254549, 545740, 252160, 232243, 235191,   268066, 249558, 199189, 211060, 227038,
               230385, 275819, 224058, 201621, 169477,   198796, 263699, 261263, 196425, 3196)
# From running awk and wc on the post-filtering vcf file
SNPdensChr <- SNPperChr/ChrSizeMb
SNPdensity <- as.data.frame(cbind(Chr, ChrSizeMb, SNPperChr, SNPdensChr))


# Total
ggplot(data = SNPdensity, aes(x = Chr, y = SNPperChr)) + 
  geom_point() +
  labs(x = "Chromosome", y = "Number of SNPs") +
  ggtitle("Number of SNPs per Chromosome") +
  theme(plot.title = element_text(hjust = 0.5)) 
ggplot(data = SNPdensity, aes(x = Chr, y = SNPdensChr/14)) + 
  geom_point() +
  labs(x = "Chromosome", y = "Number of SNPs per MB") +
  ggtitle("Number of SNPs per MB of Chromosome") +
  theme(plot.title = element_text(hjust = 0.5))
# The latter shows the SNP density, the former SNP count



###################
# Calculate F_ROH #
###################

# Total length of ROHs (short + long) / genome size
# One value per individual
c <- 1
for (i in dta$Sire){
  dta$LROH[c] <- sum(ROH$KB[ROH$IID==i])
  dta$LROHlong[c] <- sum(ROH$KB[ROH$IID==i&ROH$KB>=1000])
  dta$LROHm_l[c] <- sum(ROH$KB[ROH$IID==i&ROH$KB>=500])
  # dta$LROHstrict[c] <- sum(ROHStrict$KB[ROHStrict$IID == i])
  c <- c + 1
}

c <- 1
for (i in dads$Sire){
  dads$LROH[c] <- sum(ROH$KB[ROH$IID == i])
  dads$LROHlong[c] <- sum(ROH$KB[ROH$IID==i&ROH$KB>=1000])
  dads$LROHm_l[c] <- sum(ROH$KB[ROH$IID==i&ROH$KB>=500])
  # dads$LROHstrict[c] <- sum(ROHStrict$KB[ROHStrict$IID == i])
  c <- c + 1
}

GenSize <- 2.06807*10^6
# Genome size in Kb
# Data gotten from https://www.ncbi.nlm.nih.gov/genome/?term=Coregonus

dta$FROH <- dta$LROH/ GenSize
dads$FROH <- dads$LROH/ GenSize
dta$FROHlong <- dta$LROHlong/ GenSize
dads$FROHlong <- dads$LROHlong/ GenSize
dta$FROHm_l <- dta$LROHm_l/ GenSize
dads$FROHm_l <- dads$LROHm_l/ GenSize
# dta$FROHstrict <- dta$LROHstrict/ GenSize
# dads$FROHstrict <- dads$LROHstrict/ GenSize


#######################
# Assign ROH category #
#######################

# ROH catergorization, Ceballos et al, 2018
# Short           100Kb to 500Kb      "10s to a few 100s"
# Medium          500Kb to 1Mb  +     "a few 100s to 1-2MB"
# Long            1Mb or more         "longer than 1-2Mb"

for(i in c(1:nrow(ROH))){
  if (ROH$KB[i] >= 1000.0){
    ROH$CAT[i] <- "long"
  }else if (ROH$KB[i] >= 500.0){
    ROH$CAT[i] <- "medium"
  }else{
    ROH$CAT[i] <- "short"
  }
}
ROH$CAT <- as.factor(ROH$CAT)
summary(ROH$CAT)

# for(i in c(1:nrow(ROHStrict))){
#   if (ROHStrict$KB[i] >= 1000.0){
#     ROHStrict$CAT[i] <- "long"
#   }else if (ROHStrict$KB[i] >= 500.0){
#     ROHStrict$CAT[i] <- "medium"
#   }else{
#     ROHStrict$CAT[i] <- "short"
#   }
# }
# ROHStrict$CAT <- as.factor(ROHStrict$CAT)
# summary(ROHStrict$CAT)

summary(ROH$KB)

ggplot(data = ROH, aes(x = KB)) + 
  geom_histogram(binwidth = 100) +
  labs(x = "ROH Length", y = "Frequency") +
  ggtitle("ROH Frequency") +
  theme(plot.title = element_text(hjust = 0.5))

# ROH2 <- subset(ROH, KB >= 500)
# ggplot(data = ROH2, aes(x = KB)) + 
#   geom_histogram(binwidth = 100) +
#   labs(x = "ROH Length", y = "Frequency") +
#   ggtitle("Medium and Long ROH (> 0.5MB) Frequency") +
#   theme(plot.title = element_text(hjust = 0.5))

ROH3 <- subset(ROH, KB >= 1000)
ggplot(data = ROH3, aes(x = KB)) + 
  geom_histogram(binwidth = 100) +
  labs(x = "ROH Length", y = "Frequency") +
  ggtitle("Long ROH (> 1MB) Frequency") +
  theme(plot.title = element_text(hjust = 0.5))





#############################
# Graph ROHs per Chromosome #
#############################

ChrInd <- sort(rep(seq(1, 40, 1), 14), decreasing = F)
Sire <- rep(c(1102, 1103, 1105, 1106, 1107, 1108, 1109, 1110, 1111, 1112, 1113, 1114, 1115, 1117), 40)
ChrSize <- c(rep(93.46, 14), rep(43.33, 14), rep(42.76, 14), rep(92.22, 14), rep(45.59, 14),   rep(43.69, 14), rep(65.39, 14), rep(68.14, 14), rep(65.19, 14), rep(60.47, 14),
             rep(63.18, 14), rep(63.88, 14), rep(57.74, 14), rep(47.26, 14), rep(55.64, 14),   rep(54.03, 14), rep(54.22, 14), rep(51.95, 14), rep(59.91, 14), rep(54.34, 14),
             rep(52.95, 14), rep(56.86, 14), rep(52.02, 14), rep(50.33, 14), rep(51.03, 14),   rep(50.92, 14), rep(48.68, 14), rep(46.67, 14), rep(48.98, 14), rep(48.68, 14),
             rep(48.45, 14), rep(44.62, 14), rep(44.66, 14), rep(40.73, 14), rep(42.61, 14),   rep(42.01, 14), rep(43.66, 14), rep(36.77, 14), rep(33.96, 14), rep(01.09, 14))
ChrSize <- 1000*ChrSize

# Initialize the vectors we'll need to fill
nROH <- rep(0, length(ChrInd))
nROHlong <- rep(0, length(ChrInd))
nROHml <- rep(0, length(ChrInd))
lROH <- rep(0, length(ChrInd))
lROHlong <- rep(0, length(ChrInd))
lROHml <- rep(0, length(ChrInd))

# Make the data frame
ROHdensity <- as.data.frame(cbind(ChrInd, ChrSize, Sire, nROH, nROHlong, nROHml, lROH, lROHlong, lROHml))
ROHdensity$Sire <- as.factor(ROHdensity$Sire)

for (c in unique(ROH$CHR)){
  for (i in unique(ROH$IID)){
    # Count ROH occurence for each chromosome and individual
    ROHdensity$nROH[ROHdensity$ChrInd==c&ROHdensity$Sire==i] <- length(ROH$KB[ROH$IID==i&ROH$CHR==c])
    ROHdensity$nROHlong[ROHdensity$ChrInd==c&ROHdensity$Sire==i] <- length(ROH$KB[ROH$IID==i&ROH$CHR==c&ROH$KB>=1000])
    ROHdensity$nROHml[ROHdensity$ChrInd==c&ROHdensity$Sire==i] <- length(ROH$KB[ROH$IID==i&ROH$CHR==c&ROH$KB>=500])
    # Count ROH length over each chromosome and individual
    ROHdensity$lROH[ROHdensity$ChrInd==c&ROHdensity$Sire==i] <- sum(ROH$KB[ROH$IID==i&ROH$CHR==c])
    ROHdensity$lROHlong[ROHdensity$ChrInd==c&ROHdensity$Sire==i] <- sum(ROH$KB[ROH$IID==i&ROH$CHR==c&ROH$KB>=1000])
    ROHdensity$lROHml[ROHdensity$ChrInd==c&ROHdensity$Sire==i] <- sum(ROH$KB[ROH$IID==i&ROH$CHR==c&ROH$KB>=500])
  }
}


### Number of ROH per chromosome
ggplot(data = ROHdensity, aes(x = ChrInd, y = nROH, fill = Sire)) +
  geom_bar(position = "dodge", stat = "identity") +
  scale_fill_viridis(discrete = T) +
  labs(x = "Chromosome", y = "ROH Count") +
  ggtitle("Number of ROHs per Chromosome") +
  theme(plot.title = element_text(hjust = 0.5))

ggplot(data = ROHdensity, aes(x = ChrInd, y = nROHlong, fill = Sire)) +
  geom_bar(position = "dodge", stat = "identity") +
  scale_fill_viridis(discrete = T) +
  labs(x = "Chromosome", y = "Long ROH (> 1MB) Count") +
  ggtitle("Number of Long ROHs per Chromosome") +
  theme(plot.title = element_text(hjust = 0.5))

# ggplot(data = ROHdensity, aes(x = ChrInd, y = nROHml, fill = Sire)) +
#   geom_bar(position = "dodge", stat = "identity") +
#   scale_fill_viridis(discrete = T) +
#   labs(x = "Chromosome", y = "ROH Count") +
#   ggtitle("Number of ROHs per Chromosome") +
#   theme(plot.title = element_text(hjust = 0.5))





### Length of ROH per chromosome
ggplot(data = ROHdensity, aes(x = ChrInd, y = (lROH/ChrSize), fill = Sire)) +
  geom_bar(position = "dodge", stat = "identity") +
  scale_fill_viridis(discrete = T) +
  labs(x = "Chromosome", y = "ROH Count per MB") +
  ggtitle("Proportion of Chromosomes Covered by ROHs") +
  theme(plot.title = element_text(hjust = 0.5))

ggplot(data = ROHdensity, aes(x = ChrInd, y = (lROHlong/ChrSize), fill = Sire)) +
  geom_bar(position = "dodge", stat = "identity") +
  scale_fill_viridis(discrete = T) +
  labs(x = "Chromosome", y = "Long ROH (> 1MB) Count per MB") +
  ggtitle("Proportion of Chromosomes Covered by long ROHs") +
  theme(plot.title = element_text(hjust = 0.5))

# ggplot(data = ROHdensity, aes(x = ChrInd, y = (lROHml/ChrSize), fill = Sire)) +
#   geom_bar(position = "dodge", stat = "identity") +
#   #scale_x_continuous(breaks = seq(1, 40, 1))+
#   scale_fill_viridis(discrete = T) +
#   labs(x = "Chromosome", y = "Medium and Long ROH (> 0.5MB) Count per MB") +
#   ggtitle("Proportion of Chromosomes Covered by medium and long ROHs") +
#   theme(plot.title = element_text(hjust = 0.5))

 



#############################
# ROH density per Chromosome #
##############################

Chr <- seq(1, 40, 1)
ChrSizeMb <- c(93.46, 43.33, 42.76, 92.22, 45.59,   43.69, 65.39, 68.14, 65.19, 60.47,
               63.18, 63.88, 57.74, 47.26, 55.64,   54.03, 54.22, 51.95, 59.91, 54.34,
               52.95, 56.86, 52.02, 50.33, 51.03,   50.92, 48.68, 46.67, 48.98, 48.68,
               48.45, 44.62, 44.66, 40.73, 42.61,   42.01, 43.66, 36.77, 33.96, 01.09)

ROHdens <- matrix(nrow = 40, ncol = 1)
for(i in c(1:40)){
  ROHdens[i] <- sum(ROH$KB[ROH$CHR == i])
}

ROHbyChr <- cbind(Chr, ChrSizeMb, ROHdens)
ROHbyChr <- as.data.frame(ROHbyChr)
colnames(ROHbyChr) <- c("Chr", "ChrSizeMb", "ROHdens")

for (i in ROH$CHR) {
  nlong <- nrow(ROH[ROH$CHR == i & ROH$CAT == "long",])
  nmed <- nrow(ROH[ROH$CHR == i & ROH$CAT == "medium",])
  nshort <- nrow(ROH[ROH$CHR == i & ROH$CAT == "short",])
  ROHbyChr$nlong[i] <- nlong
  ROHbyChr$nmedium[i] <- nmed
  ROHbyChr$nshort[i] <- nshort
}

ROHbyChr$Total <- ROHbyChr$nlong + ROHbyChr$nmedium + ROHbyChr$nshort



##########################################
### ROH categories for each individual ###
##########################################

IndID <- unique(ROH$IID)
FROHid <- rep(0, length(IndID))
nlong <- rep(0, length(IndID))
nmedium <- rep(0, length(IndID))
nshort <- rep(0, length(IndID))
ntotal <- rep(0, length(IndID))
ROHperInd <- as.data.frame(cbind(IndID, FROHid, ntotal, nlong, nmedium, nshort))

for (i in IndID){
  ROHperInd$FROHid[ROHperInd$IndID==i] <- dads$FROHlong[dads$Sire==i]
  ROHperInd$nlong[ROHperInd$IndID==i] <- length(ROH$SNP1[ROH$CAT=="long"&ROH$IID==i])
  ROHperInd$nmedium[ROHperInd$IndID==i] <- length(ROH$SNP1[ROH$CAT=="medium"&ROH$IID==i])
  ROHperInd$nshort[ROHperInd$IndID==i] <- length(ROH$SNP1[ROH$CAT=="short"&ROH$IID==i])
  ROHperInd$ntotal[ROHperInd$IndID==i] <- length(ROH$SNP1[ROH$IID==i])
}





#################################
# Calculate Nucleotide Dversity #
#################################

# Calculate nucleotide diversity (pi)
# Nei and Li, 1973
# pi = sum_i^j(x_i x_j pi_{ij})


# Unsure how to best do this 



###############################################
# Function for Generalized Mixed Effect Model #
###############################################

makeGLMModelSireDam <- function(dta, respVar, fixedPar1, fixedPar2){
  # Check which variables were filled in
  if(missing(fixedPar2)){
    if(missing(fixedPar1)){
      if(missing(respVar)){
        if(missing(dta)){
          stop('You need to include data, a response variable and one or two fixed explanatory variables')
        }else{
          stop('You need to include a response variable and one or two fixed explanatory variables')
        }
      }else{
        stop('You need to include one or two fixed explanatory variable')
      }
    }else{
      # Make reference model with Response Variable, Fixed Explanatory Variable, Sire and Dam
      ref <- glmer(respVar ~ fixedPar1 + (1|Sire) + (1|Dam), data=dta, family=binomial, control=glmerControl(check.nlev.gtr.1="ignore"), na.action = na.exclude) # reference model
      # Make alternate models each excluding one of the explanatory variables, no interaction models
      mod1 <- glmer(respVar ~ fixedPar1 + (1|Sire), data=dta, family=binomial, control=glmerControl(check.nlev.gtr.1="ignore"), na.action = na.exclude) # dam
      mod2 <- glmer(respVar ~ fixedPar1 + (1|Dam), data=dta, family=binomial, control=glmerControl(check.nlev.gtr.1="ignore"), na.action = na.exclude) # sire
      mod3 <- glmer(respVar ~ (1|Sire) + (1|Dam), data=dta, family=binomial, control=glmerControl(check.nlev.gtr.1="ignore"), na.action = na.exclude) # fixedPar1
      # Likelihood ratio test
      # Assesses the goodness of fit of two competing statistical models based on the ratio of their likelihoods
      an1 <- anova (ref, mod1) # dam
      an2 <- anova (ref, mod2) # sire
      an3 <- anova (ref, mod3) # fixedPar1
      # Names for output
      name1 <- paste(fixedPar1, "Effect")
      name2 <- paste("Likel", fixedPar1, "Ref")
      # Create output with all the models and likelihood ratio
      modelResults <- list("Reference Model" = ref,
                           "Dam Effect" = mod1,
                           "Sire Effect" = mod2,
                           name1 = mod3,
                           "Likel Dam Ref" = an1,
                           "Likel Sire Ref" = an2,
                           name2 = an3)
    }
  }else{
    # Make reference model with Response Variable, Fixed Explanatory Variable, Sire and Dam
    ref <-glmer(respVar~fixedPar1 + fixedPar2 + (1|Sire) + (1|Dam), data=dta, family=binomial, control=glmerControl(check.nlev.gtr.1="ignore"), na.action = na.exclude) # reference model
    # Make alternate models each excluding one of the explanatory variables, no interaction models
    mod1 <-glmer(respVar ~ fixedPar1 + fixedPar2 + (1|Sire), data=dta, family=binomial, control=glmerControl(check.nlev.gtr.1="ignore"), na.action = na.exclude) # Dam
    mod2 <-glmer(respVar ~ fixedPar1 + fixedPar2 + (1|Dam), data=dta, family=binomial, control=glmerControl(check.nlev.gtr.1="ignore"), na.action = na.exclude) # Sire
    mod3 <-glmer(respVar ~ fixedPar2 + (1|Sire) + (1|Dam), data=dta, family=binomial, control=glmerControl(check.nlev.gtr.1="ignore"), na.action = na.exclude) # fixedPar1
    mod4 <-glmer(respVar ~ fixedPar1 + (1|Sire) + (1|Dam), data=dta, family=binomial, control=glmerControl(check.nlev.gtr.1="ignore"), na.action = na.exclude) # fixedPar2
    # Make alternate models that take interactions into account
    mod5 <-glmer(respVar ~ fixedPar1 + fixedPar2 + (fixedPar1|Sire) + (1|Dam), data=dta, family=binomial, control=glmerControl(check.nlev.gtr.1="ignore"), na.action = na.exclude) # fixedPar1 + Sire
    mod6 <-glmer(respVar ~ fixedPar1 + fixedPar2 + fixedPar1:fixedPar2 + (1|Sire) + (1|Dam), data=dta, family=binomial, control=glmerControl(check.nlev.gtr.1="ignore"), na.action = na.exclude) # fixedPar1 + fixedPar2
    # Likelihood ratio test
    # Assesses the goodness of fit of two competing statistical models based on the ratio of their likelihoods
    an1 <- anova (ref, mod1) # Dam
    an2 <- anova (ref, mod2) # Sire
    an3 <- anova (ref, mod3) # fixedPar1
    an4 <- anova (ref, mod4) # fixedPar2
    an5 <- anova (ref, mod5) # fixedPar1 x Sire
    an6 <- anova (ref, mod6) # fixedPar1 x fixedPar2
    # Names for output
    name1 <- paste(fixedPar1, "Effect")
    name2 <- paste(fixedPar2, "Effect")
    name3 <- paste(fixedPar1, "x Sire Effect")
    name4 <- paste(fixedPar1, "x", fixedPar2, "Effect")
    name5 <- paste("Likel", fixedPar1, "Ref")
    name6 <- paste("Likel", fixedPar2, "Ref")
    name7 <- paste("Likel", fixedPar1, "x Sire Ref")
    name8 <- paste("Likel", fixedPar1, "x", fixedPar2, "Ref")
    # Create output with all the models and likelihood ratio
    modelResults <- list("Reference Model" = ref,
                         "Dam Effect" = mod1,
                         "Sire Effect" = mod2,
                         name1 = mod3,
                         name2 = mod4,
                         name3 = mod5,
                         name4 = mod6,
                         "Likel Dam Ref" = an1,
                         "Likel Sire Ref" = an2,
                         name5 = an3,
                         name6 = an4,
                         name7 = an5,
                         name8 = an6)
    # modelResults <- list("Reference Model" = ref,
    #                      "Dam Effect" = mod1,
    #                      "Sire Effect" = mod2,
    #                      name1 = mod3,
    #                      name2 = mod4,
    #                      name3 = mod5,
    #                      name4 = mod6)
    # names(modelResults) <- c("Reference Model", "Dam Effect", "Sire Effect", 
    #                          name1, name2, name3, name4, "Likel Dam Ref",
    #                          "Likel Sire Ref", name5, name6, name7, name8)
  }
  return(modelResults)
}



###################################
# Function for Mixed Effect Model #
###################################

makeLMModelSireDam <- function(dta, respVar, fixedPar1, fixedPar2){
  # Check which variables were filled in
  if(missing(fixedPar2)){
    if(missing(fixedPar1)){
      if(missing(respVar)){
        if(missing(dta)){
          stop('You need to include data, a response variable and one or two fixed explanatory variables')
        }else{
          stop('You need to include a response variable and one or two fixed explanatory variables')
        }
      }else{
        stop('You need to include one or two fixed explanatory variable')
      }
    }else{
      # Make reference model with Response Variable, Fixed Explanatory Variable, Sire and Dam
      ref <- lmer(respVar~fixedPar1 + (1|Sire) + (1|Dam), data= dta, REML = FALSE, control=lmerControl(check.nlev.gtr.1="ignore")) # reference model
      # Make alternate models each excluding one of the explanatory variables, no interaction models
      mod1 <- lmer(respVar ~ fixedPar1 + (1|Sire), data = dta, REML = FALSE, control=lmerControl(check.nlev.gtr.1="ignore")) # dam
      mod2 <- lmer(respVar ~ fixedPar1 + (1|Dam), data = dta, REML = FALSE, control=lmerControl(check.nlev.gtr.1="ignore")) # sire
      mod3 <- lmer(respVar ~ (1|Sire) + (1|Dam), data = dta, REML = FALSE, control=lmerControl(check.nlev.gtr.1="ignore")) # treatment
      # Likelihood ratio test
      # Assesses the goodness of fit of two competing statistical models based on the ratio of their likelihoods
      an1 <- anova (ref, mod1) # dam
      an2 <- anova (ref, mod2) # sire
      an3 <- anova (ref, mod3) # fixedPar1
      # Names for output
      name1 <- paste(fixedPar1, "Effect")
      name2 <- paste("Likel", fixedPar1, "Ref")
      # Create output with all the models and likelihood ratio
      modelResults <- list("Reference Model" = ref,
                           "Dam Effect" = mod1,
                           "Sire Effect" = mod2,
                           name1 = mod3,
                           "Likel Dam Ref" = an1,
                           "Likel Sire Ref" = an2,
                           name2 = an3)
    }
  }else{
    # Make reference model with Response Variable, Fixed Explanatory Variable, Sire and Dam
    ref <- lmer(respVar~fixedPar1 + fixedPar2 + (1|Sire) + (1|Dam), data= dta, REML = FALSE, control=lmerControl(check.nlev.gtr.1="ignore")) # reference model
    # Make alternate models each excluding one of the explanatory variables, no interaction models
    mod1 <- lmer(respVar~fixedPar1 + fixedPar2 + (1|Sire), data= dta, REML = FALSE, control=lmerControl(check.nlev.gtr.1="ignore")) # dam
    mod2 <- lmer(respVar~fixedPar1 + fixedPar2 + (1|Dam), data= dta, REML = FALSE, control=lmerControl(check.nlev.gtr.1="ignore")) # sire
    mod3 <- lmer(respVar~fixedPar2 + (1|Sire) + (1|Dam), data= dta, REML = FALSE, control=lmerControl(check.nlev.gtr.1="ignore")) # treatment
    mod4 <- lmer(respVar~fixedPar1 + (1|Sire) + (1|Dam), data= dta, REML = FALSE, control=lmerControl(check.nlev.gtr.1="ignore")) # fixedPar2
    # Make alternate models that take interactions into account
    mod5 <- lmer(respVar~fixedPar1 + fixedPar2 + (fixedPar1|Sire) + (1|Dam), data= dta, REML = FALSE, control=lmerControl(check.nlev.gtr.1="ignore")) # treatment-sire
    mod6 <- lmer(respVar~fixedPar1 + fixedPar2 + fixedPar1:fixedPar2 + (1|Sire) + (1|Dam), data= dta, REML = FALSE, control=lmerControl(check.nlev.gtr.1="ignore")) # treatment fixedPar2
    # Likelihood ratio test
    # Assesses the goodness of fit of two competing statistical models based on the ratio of their likelihoods
    an1 <- anova (ref, mod1) # Dam
    an2 <- anova (ref, mod2) # Sire
    an3 <- anova (ref, mod3) # fixedPar1
    an4 <- anova (ref, mod4) # fixedPar2
    an5 <- anova (ref, mod5) # fixedPar1 x Sire
    an6 <- anova (ref, mod6) # fixedPar1 x fixedPar2
    # Names for output
    name1 <- paste(fixedPar1, "Effect")
    name2 <- paste(fixedPar2, "Effect")
    name3 <- paste(fixedPar1, "x Sire Effect")
    name4 <- paste(fixedPar1, "x", fixedPar2, "Effect")
    name5 <- paste("Likel", fixedPar1, "Ref")
    name6 <- paste("Likel", fixedPar2, "Ref")
    name7 <- paste("Likel", fixedPar1, "x Sire Ref")
    name8 <- paste("Likel", fixedPar1, "x", fixedPar2, "Ref")
    # Create output with all the models and likelihood ratio
    modelResults <- list("Reference Model" = ref,
                         "Dam Effect" = mod1,
                         "Sire Effect" = mod2,
                         name1 = mod3,
                         name2 = mod4,
                         name3 = mod5,
                         name4 = mod6,
                         "Likel Dam Ref" = an1,
                         "Likel Sire Ref" = an2,
                         name5 = an3,
                         name6 = an4,
                         name7 = an5,
                         name8 = an6)
  }
  return(modelResults)
}


########################
### MODELS AND PLOTS ###
########################


#################
# Outliers test #
#################

### Eggs volume before hatching ###
plot(dta$Embryo_volume , type="b")
out <- outlier(dta$Embryo_volume) # returns the most extreme value:
# largest diff. between this point and the mean of the data
grubbs.test(dta$Embryo_volume)### if it is an outlier the test is significant!


dta$Embryo_volume[dta$Embryo_volume == 8.359220854]<- NA
dta$Embryo_volume[dta$Embryo_volume == 8.309522882]<- NA


### Length at hatching ###
plot(dta$Length_D0 , type="b")
out <- outlier(dta$Length_D0) # returns the most extreme value:
# largest diff. between this point and the mean of the data
grubbs.test(dta$Length_D0)### if it is an outlier the test is significant!

dta$Length_D0[dta$Length_D0 == 5.408] <- NA
dta$Length_D0[dta$Length_D0 == 5.724] <- NA
dta$Length_D0[dta$Length_D0 == 5.833] <- NA
dta$Length_D0[dta$Length_D0 == 5.954] <- NA
dta$Length_D0[dta$Length_D0 == 6.013] <- NA
dta$Length_D0[dta$Length_D0 == 6.138] <- NA
dta$Length_D0[dta$Length_D0 == 6.225] <- NA
dta$Length_D0[dta$Length_D0 == 6.249] <- NA
dta$Length_D0[dta$Length_D0 == 6.299] <- NA
dta$Length_D0[dta$Length_D0 == 6.3] <- NA
dta$Length_D0[dta$Length_D0 == 6.325] <- NA
dta$Length_D0[dta$Length_D0 == 6.456] <- NA
dta$Length_D0[dta$Length_D0 == 6.481] <- NA
dta$Length_D0[dta$Length_D0 == 6.52] <- NA
dta$Length_D0[dta$Length_D0 == 6.662] <- NA
dta$Length_D0[dta$Length_D0 == 6.825] <- NA
dta$Length_D0[dta$Length_D0 == 6.885] <- NA
dta$Length_D0[dta$Length_D0 == 6.225] <- NA
dta$Length_D0[dta$Length_D0 == 6.881] <- NA
dta$Length_D0[dta$Length_D0 == 6.932] <- NA
dta$Length_D0[dta$Length_D0 == 6.954] <- NA
dta$Length_D0[dta$Length_D0 == 6.959] <- NA
dta$Length_D0[dta$Length_D0 == 6.964] <- NA
dta$Length_D0[dta$Length_D0 == 6.98] <- NA
dta$Length_D0[dta$Length_D0 == 7.075] <- NA
dta$Length_D0[dta$Length_D0 == 7.111] <- NA
dta$Length_D0[dta$Length_D0 == 7.156] <- NA
dta$Length_D0[dta$Length_D0 == 7.167] <- NA
dta$Length_D0[dta$Length_D0 == 7.183] <- NA
dta$Length_D0[dta$Length_D0 == 7.189] <- NA


### Length day21 ###
plot(dta$Length_D21 , type="b")
out <- outlier(dta$Length_D21) # returns the most extreme value:
# largest diff. between this point and the mean of the data
grubbs.test(dta$Length_D21)### if it is an outlier the test is significant!

dta$Length_D21[dta$Length_D21 == 6.257] <- NA
dta$Length_D21[dta$Length_D21 == 6.505] <- NA
dta$Length_D21[dta$Length_D21 == 6.509] <- NA
dta$Length_D21[dta$Length_D21 == 7.239] <- NA
dta$Length_D21[dta$Length_D21 == 7.251] <- NA
dta$Length_D21[dta$Length_D21 == 7.286] <- NA
dta$Length_D21[dta$Length_D21 == 7.348] <- NA
dta$Length_D21[dta$Length_D21 == 7.473] <- NA
dta$Length_D21[dta$Length_D21 == 7.778] <- NA
dta$Length_D21[dta$Length_D21 == 8] <- NA
dta$Length_D21[dta$Length_D21 == 8.007] <- NA
dta$Length_D21[dta$Length_D21 == 8.031] <- NA
dta$Length_D21[dta$Length_D21 == 8.051] <- NA


### Growth rate D21 ###
plot(dta$Growth_rateD21 , type="b")
out <- outlier(dta$Growth_rateD21) # returns the most extreme value:
# largest diff. between this point and the mean of the data
grubbs.test(dta$Growth_rateD21)### if it is an outlier the test is significant!

dta$Growth_rateD21[dta$Growth_rateD21 > 0.89] <- NA
dta$Growth_rateD21[dta$Growth_rateD21 > 0.76] <- NA
dta$Growth_rateD21[dta$Growth_rateD21 > 0.73] <- NA
dta$Growth_rateD21[dta$Growth_rateD21 > 0.68] <- NA
dta$Growth_rateD21[dta$Growth_rateD21 > 0.60] <- NA
dta$Growth_rateD21[dta$Growth_rateD21 > 0.58] <- NA
dta$Growth_rateD21[dta$Growth_rateD21 > 0.56] <- NA
dta$Growth_rateD21[dta$Growth_rateD21 > 0.55] <- NA
dta$Growth_rateD21[dta$Growth_rateD21 > 0.50] <- NA
dta$Growth_rateD21[dta$Growth_rateD21 > 0.46] <- NA
dta$Growth_rateD21[dta$Growth_rateD21 > 0.45] <- NA
dta$Growth_rateD21[dta$Growth_rateD21 > 0.44] <- NA
dta$Growth_rateD21[dta$Growth_rateD21 > 0.43] <- NA
dta$Growth_rateD21[dta$Growth_rateD21 > 0.42] <- NA
dta$Growth_rateD21[dta$Growth_rateD21 > 0.41] <- NA
dta$Growth_rateD21[dta$Growth_rateD21 > 0.40] <- NA
dta$Growth_rateD21[dta$Growth_rateD21 > 0.39] <- NA
dta$Growth_rateD21[dta$Growth_rateD21 > 0.38] <- NA
dta$Growth_rateD21[dta$Growth_rateD21 > 0.37] <- NA



#############
# Offspring #
#############


### Yolk sac volume and ROH ###
# makeLMModelSireDam(dta, dta$YS_volume, dta$Treatment_1, dta$FROHm_l)
makeLMModelSireDam(dta, dta$YS_volume, dta$Treatment_1, dta$FROHlong)

### Hatching Time and ROH ###
# makeLMModelSireDam(dta, dta$Hatching_days, dta$Treatment_1, dta$FROHm_l)
makeLMModelSireDam(dta, dta$Hatching_days, dta$Treatment_1, dta$FROHlong)

### Egg Volume before hatching and ROH ###
# makeLMModelSireDam(dta, dta$Embryo_volume, dta$Treatment_1, dta$FROHm_l)
makeLMModelSireDam(dta, dta$Embryo_volume, dta$Treatment_1, dta$FROHlong)

### Length at Hatching and ROH ###
# makeLMModelSireDam(dta, dta$Length_D0, dta$Treatment_1, dta$FROHm_l)
makeLMModelSireDam(dta, dta$Length_D0, dta$Treatment_1, dta$FROHlong)

### Length on Final Day and ROH ###
# makeLMModelSireDam(dta, dta$Length_D21, dta$Treatment_1, dta$FROHm_l)
makeLMModelSireDam(dta, dta$Length_D21, dta$Treatment_1, dta$FROHlong)

### Growth Rate at day 21 and ROH ###
# makeLMModelSireDam(dta, dta$Growth_rateD21, dta$Treatment_1, dta$FROHm_l)
makeLMModelSireDam(dta, dta$Growth_rateD21, dta$Treatment_1, dta$FROHlong)






### Treatment_1 effect on egg volume ###
list<- list(dta$Sibgroup, dta$Treatment_1)
MeansOfMeansEmbryo <- tapply(dta$Embryo_volume, list, mean, na.rm=TRUE)
MeanCEmbryo<- mean(MeansOfMeansEmbryo[,1], na.rm=TRUE) #Mean Control
MeanPFEmbryo<- mean(MeansOfMeansEmbryo[,2], na.rm=TRUE) #Mean PF
vector.meanEmbryo <- c(MeanCEmbryo, MeanPFEmbryo)
sd.C.Embryo <- sd(dta$Embryo_volume[dta$Treatment_1=="C"], na.rm=TRUE) #SD Control
sd.PF.Embryo <- sd(dta$Embryo_volume[dta$Treatment_1=="PF2"], na.rm=TRUE)# SD PF
conf.int.C.Embryo <- qnorm(0.975)*sd.C.Embryo/sqrt(length(dta$Embryo_volume[dta$Treatment_1=="C"]))
conf.int.PF.Embryo <- qnorm(0.975)*sd.PF.Embryo /sqrt(length(dta$Embryo_volume[dta$Treatment_1=="PF2"]))
vector.conf.int <- c(conf.int.C.Embryo,conf.int.PF.Embryo ) 
barplot2(vector.meanEmbryo,plot.ci=T, ylim= c(3, 5), names=c("C", "PF2"), ci.u=vector.meanEmbryo+vector.conf.int, ci.l=vector.meanEmbryo-vector.conf.int, ylab="Eggs volume", xpd=FALSE)

### Treatment_1 effect on hatching time ###
list<- list(dta$Sibgroup, dta$Treatment_1)
MeansOfMeansHT <- tapply(dta$Hatching_days, list, mean, na.rm=TRUE)
MeanCHT<- mean(MeansOfMeansHT[,1], na.rm=TRUE) #Mean HT Control
MeanPFHT<- mean(MeansOfMeansHT[,2], na.rm=TRUE) #Mean HT PF
vector.meanHT <- c(MeanCHT, MeanPFHT)
sd.C.HT <- sd(dta$Hatching_days[dta$Treatment_1=="C"], na.rm=TRUE) #SD Control
sd.PF.HT <- sd(dta$Hatching_days[dta$Treatment_1=="PF2"], na.rm=TRUE)# SD PF
conf.int.C.HT <- qnorm(0.975)*sd.C.HT/sqrt(length(dta$Hatching_days[dta$Treatment_1=="C"]))
conf.int.PF.HT <- qnorm(0.975)*sd.PF.HT /sqrt(length(dta$Hatching_days[dta$Treatment_1=="PF2"]))
vector.conf.int <- c(conf.int.C.HT,conf.int.PF.HT ) 
barplot2(vector.meanHT,plot.ci=T, ylim= c(85, 94), names=c("C", "PF2"), ci.u=vector.meanHT+vector.conf.int, ci.l=vector.meanHT-vector.conf.int, ylab="Hatching time (days)", xpd=FALSE)

### Treatment_1 effect on Yolk sac volume at hatching ###
list<- list(dta$Sibgroup, dta$Treatment_1)
MeansOfMeansYS <- tapply(dta$YS_volume, list, mean, na.rm=TRUE)
MeanCYS<- mean(MeansOfMeansYS[,1], na.rm=TRUE) #Mean YS Control
MeanPFYS<- mean(MeansOfMeansYS[,2], na.rm=TRUE) #Mean YS PF
vector.meanYS <- c(MeanCYS, MeanPFYS)
sd.C.YS <- sd(dta$YS_volume[dta$Treatment_1=="C"], na.rm=TRUE) #SD Control
sd.PF.YS <- sd(dta$YS_volume[dta$Treatment_1=="PF2"], na.rm=TRUE)# SD PF
conf.int.C.YS <- qnorm(0.975)*sd.C.YS/sqrt(length(dta$YS_volume[dta$Treatment_1=="C"]))
conf.int.PF.YS <- qnorm(0.975)*sd.PF.YS /sqrt(length(dta$YS_volume[dta$Treatment_1=="PF2"]))
vector.conf.int <- c(conf.int.C.YS,conf.int.PF.YS ) 
barplot2(vector.meanYS,plot.ci=T, ylim= c(0.5, 0.7), names=c("C", "PF2"), ci.u=vector.meanYS+vector.conf.int, ci.l=vector.meanYS-vector.conf.int, ylab="Yolk sac volume", xpd=FALSE)

### Treatment_1 effect on Length at hatching###
list<- list(dta$Sibgroup, dta$Treatment_1)
MeansOfMeansLD0 <- tapply(dta$Length_D0, list, mean, na.rm=TRUE)
MeanCLD0<- mean(MeansOfMeansLD0[,1], na.rm=TRUE) #Mean LD0 Control
MeanPFLD0<- mean(MeansOfMeansLD0[,2], na.rm=TRUE) #Mean LD0 PF
vector.meanLD0 <- c(MeanCLD0, MeanPFLD0)
sd.C.LD0 <- sd(dta$Length_D0[dta$Treatment_1=="C"], na.rm=TRUE) #SD Control
sd.PF.LD0 <- sd(dta$Length_D0[dta$Treatment_1=="PF2"], na.rm=TRUE)# SD PF
conf.int.C.LD0 <- qnorm(0.975)*sd.C.LD0/sqrt(length(dta$Length_D0[dta$Treatment_1=="C"]))
conf.int.PF.LD0 <- qnorm(0.975)*sd.PF.LD0 /sqrt(length(dta$Length_D0[dta$Treatment_1=="PF2"]))
vector.conf.int <- c(conf.int.C.LD0,conf.int.PF.LD0 ) 
barplot2(vector.meanLD0,plot.ci=T, ylim= c(9, 9.5), names=c("C", "PF2"), ci.u=vector.meanLD0+vector.conf.int, ci.l=vector.meanLD0-vector.conf.int, ylab="Length at hatching", xpd=FALSE)


### Treatment_1 effect on final length ###
list<- list(dta$Sibgroup, dta$Treatment_1)
MeansOfMeansLD21 <- tapply(dta$Length_D21, list, mean, na.rm=TRUE)
MeanCLD21<- mean(MeansOfMeansLD21[,1], na.rm=TRUE) #Mean LD21 Control
MeanPFLD21<- mean(MeansOfMeansLD21[,2], na.rm=TRUE) #Mean LD21 PF
vector.meanLD21 <- c(MeanCLD21, MeanPFLD21)
sd.C.LD21 <- sd(dta$Length_D21[dta$Treatment_1=="C"], na.rm=TRUE) #SD Control
sd.PF.LD21 <- sd(dta$Length_D21[dta$Treatment_1=="PF2"], na.rm=TRUE)# SD PF
conf.int.C.LD21 <- qnorm(0.975)*sd.C.LD21/sqrt(length(dta$Length_D21[dta$Treatment_1=="C"]))
conf.int.PF.LD21 <- qnorm(0.975)*sd.PF.LD21 /sqrt(length(dta$Length_D21[dta$Treatment_1=="PF2"]))
vector.conf.int <- c(conf.int.C.LD21,conf.int.PF.LD21 ) 
barplot2(vector.meanLD21,plot.ci=T, ylim= c(10, 10.5), names=c("C", "PF2"), ci.u=vector.meanLD21+vector.conf.int, ci.l=vector.meanLD21-vector.conf.int, ylab="Length day 21", xpd=FALSE)

### Treatment_1 effect on growth rate after 21 days ####
list<- list(dta$Sibgroup, dta$Treatment_1)
MeansOfMeansGRD21 <- tapply(dta$Growth_rateD21, list, mean, na.rm=TRUE)
MeanCGRD21<- mean(MeansOfMeansGRD21[,1], na.rm=TRUE) #Mean GRD21 Control
MeanPFGRD21<- mean(MeansOfMeansGRD21[,2], na.rm=TRUE) #Mean GRD21 PF
vector.meanGRD21 <- c(MeanCGRD21, MeanPFGRD21)
sd.C.GRD21 <- sd(dta$Growth_rateD21[dta$Treatment_1=="C"], na.rm=TRUE) #SD Control
sd.PF.GRD21 <- sd(dta$Growth_rateD21[dta$Treatment_1=="PF2"], na.rm=TRUE)# SD PF
conf.int.C.GRD21 <- qnorm(0.975)*sd.C.GRD21/sqrt(length(dta$Growth_rateD21[dta$Treatment_1=="C"]))
conf.int.PF.GRD21 <- qnorm(0.975)*sd.PF.GRD21 /sqrt(length(dta$Growth_rateD21[dta$Treatment_1=="PF2"]))
vector.conf.int <- c(conf.int.C.GRD21,conf.int.PF.GRD21 ) 
barplot2(vector.meanGRD21,plot.ci=T, ylim= c(0, 0.15), names=c("C", "PF2"), ci.u=vector.meanGRD21+vector.conf.int, ci.l=vector.meanGRD21-vector.conf.int, ylab="Growth rate day 21", xpd=FALSE)










#########
# Sires #
#########



lm_eqn <- function(x,y){
  m <- lm(y ~ x);
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2,
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3)))
  as.character(as.expression(eq));
}



#######################
# Breeding tubercules #
#######################

### Volume & long ROHs
ggplot(data = dads, aes(x = Volume, y = FROHlong)) +
  geom_point() + 
  geom_smooth(method = "lm", se=F, color = "black", formula = y ~ x) +
  labs(x = "Breeding Tubercule Volume", y = "Froh, ROH > 1Mb") +
  annotate("text", x = 0.20, y = 0.02, label = "p-value = 0.99423")

lm6 <- lm(dads$FROHlong~dads$Volume) #Control
summary(lm6)

# ### Volume & medium-long ROHs
# ggplot(data = dads, aes(x = Volume, y = FROHm_l)) +
#   geom_point() +
#   geom_smooth(method = "lm", se=F, color = "black", formula = y ~ x) +
#   labs(x = "Breeding Tubercule Volume", y = "Froh, ROH > 0.5Mb") +
#   annotate("text", x = 0.20, y = 0.06, label = "p-value = 0.9787")
# 
# lm6 <- lm(dads$FROHm_l~dads$Volume) #Control
# summary(lm6)


### Max depth & long ROHs
ggplot(data = dads, aes(x = Max_depth, y = FROHlong)) +
  geom_point() + 
  geom_smooth(method = "lm", se=F, color = "black", formula = y ~ x) +
  labs(x = "Breeding Tubercule Maximum Depth", y = "Froh, ROH > 1Mb")  +
  annotate("text", x = 0.30, y = 0.02, label = "p-value = 0.9624")

lm6 <- lm(dads$FROHlong~dads$Max_depth) #Control
summary(lm6)

# ### Max depth & medium-long ROHs
# ggplot(data = dads, aes(x = Max_depth, y = FROHm_l)) +
#   geom_point() +
#   geom_smooth(method = "lm", se=F, color = "black", formula = y ~ x) +
#   labs(x = "Breeding Tubercule Maximum Depth", y = "Froh, ROH > 0.5Mb") +
#   ggtitle("Secondary Sexual Ornaments and Froh") +
#   theme(plot.title = element_text(hjust = 0.5)) +
#   annotate("text", x = 0.30, y = 0.06, label = "p-value = 0.9562")
# 
# lm6 <- lm(dads$FROHm_l~dads$Max_depth) #Control
# summary(lm6)



#############
# Body Size #
#############


### Length & long ROHs
ggplot(data = dads, aes(x = Standard_length, y = FROHlong)) +
  geom_point() + 
  geom_smooth(method = "lm", se=F, color = "black", formula = y ~ x) +
  labs(x = "Standard Length", y = "Froh, ROH > 1Mb") +
  annotate("text", x = 250, y = 0.02, label = "p-value = 0.0536")

lm6 <- lm(dads$FROHlong~dads$Standard_length) #Control
summary(lm6)

### Length & medium-long ROHs
# ggplot(data = dads, aes(x = Standard_length, y = FROHm_l)) +
#   geom_point() +
#   geom_smooth(method = "lm", se=F, color = "black", formula = y ~ x) +
#   labs(x = "Standard Length", y = "Froh, ROH > 0.5Mb") +
#   ggtitle("Body Characteristics and Froh") +
#   theme(plot.title = element_text(hjust = 0.5)) +
#   annotate("text", x = 250, y = 0.06, label = "p-value = 0.05477")
#   
# lm6 <- lm(dads$FROHm_l~dads$Standard_length) #Control
# summary(lm6)


### Weight & long ROHs
ggplot(data = dads, aes(x = Weight, y = FROHlong)) +
  geom_point() + 
  geom_smooth(method = "lm", se=F, color = "black", formula = y ~ x) +
  labs(x = "Weight", y = "Froh, ROH > 1Mb") +
  annotate("text", x = 200, y = 0.02, label = "p-value = 0.19745")

lm6 <- lm(dads$FROHlong~dads$Weight) #Control
summary(lm6)

### Weight & medium-long ROHs
# ggplot(data = dads, aes(x = Weight, y = FROHm_l)) +
#   geom_point() +
#   geom_smooth(method = "lm", se=F, color = "black", formula = y ~ x) +
#   labs(x = "Weight", y = "Froh, ROH > 0.5Mb") +
#   ggtitle("Body Characteristics and Froh") +
#   theme(plot.title = element_text(hjust = 0.5)) +
#   annotate("text", x = 200, y = 0.06, label = "p-value = 0.2066")
# 
# lm6 <- lm(dads$FROHm_l~dads$Weight) #Control
# summary(lm6)


#########
# Sperm #
#########

### Immotiles & long ROHs
ggplot(data = dads, aes(x = Immotiles, y = FROHlong)) +
  geom_point() + 
  geom_smooth(method = "lm", se=F, color = "black", formula = y ~ x) +
  labs(x = "Immotiles", y = "Froh, ROH > 1Mb") +
  annotate("text", x = 25, y = 0.02, label = "p-value = 0.254")

lm6 <- lm(dads$FROHlong~dads$Immotiles) #Control
summary(lm6)

### Immotiles & medium-long ROHs
# ggplot(data = dads, aes(x = Immotiles, y = FROHm_l)) +
#   geom_point() +
#   geom_smooth(method = "lm", se=F, color = "black", formula = y ~ x) +
#   labs(x = "Immotiles", y = "Froh, ROH > 0.5Mb") +
#   ggtitle("Sperm Characteristics and Froh") +
#   theme(plot.title = element_text(hjust = 0.5)) +
#   annotate("text", x = 25, y = 0.06, label = "p-value = 0.2629")
# 
# lm6 <- lm(dads$FROHm_l~dads$Immotiles) #Control
# summary(lm6)


### Speed & long ROHs
ggplot(data = dads, aes(x = Speed, y = FROHlong)) +
  geom_point() + 
  geom_smooth(method = "lm", se=F, color = "black", formula = y ~ x) +
  labs(x = "Speed", y = "Froh, ROH > 1Mb") +
  annotate("text", x = 105, y = 0.02, label = "p-value = 0.354")

lm6 <- lm(dads$FROHlong~dads$Speed) #Control
summary(lm6)

### Speed & medium-long ROHs
# ggplot(data = dads, aes(x = Speed, y = FROHm_l)) +
#   geom_point() +
#   geom_smooth(method = "lm", se=F, color = "black", formula = y ~ x) +
#   labs(x = "Speed", y = "Froh, ROH > 0.5Mb") +
#   ggtitle("Sperm Characteristics and Froh") +
#   theme(plot.title = element_text(hjust = 0.5)) +
#   annotate("text", x = 105, y = 0.06, label = "p-value = 0.2361")
# 
# lm6 <- lm(dads$FROHm_l~dads$Speed) #Control
# summary(lm6)


