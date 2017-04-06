# Karl's class
#R flank vcf
#set WD

# Install Packages
install.packages("devtools")
source("http://bioconductor.org/biocLite.R")
biocLite("qvalue")
install_github("whitlock/OutFLANK")
install.packages("adgenet")

#Load the packages
library(devtools)
library(OutFLANK)
library(vcfR)
library(adegenet)

# Read in your .geno file. OutFlank requires it to be transposed (rotate and flip), so we'll do that next.
ssw.geno_in <- read.fwf("SSW_all_biallelic.MAF0.02.Miss0.8.recode.vcf.geno" , width=rep(1,24))
ssw.geno <- t(ssw.geno_in)

#check if it transposed
dim(ssw.geno)
dim(ssw.geno_in)

# Read in the meta data
ssw_meta <- read.table("ssw_healthloc.txt", header=TRUE) # read in the data
ssw_meta <- ssw_meta[order(ssw_meta$Individual),] # reorder the meta_data by Ind number
ssw_meta$Trajectory[which(ssw_meta$Trajectory == 'MM')] = NA # selects them and sets them as MM to remove the MM's from the analysis

# Now we use Out FLANK
# INPUT SNP matrix, popnames= list of the names for all ind
OF_SNPs <- MakeDiploidFSTMat(ssw.geno, locusNames = seq(1,5317,1), popNames = ssw_meta$Trajectory)
head(OF_SNPs)
dim(OF_SNPs)
# Hmin = Heterozygosity min, Qthreshold= false discovery rate
OF_out <- OutFLANK(FstDataFrame = OF_SNPs, LeftTrimFraction = 0.05, RightTrimFraction = 0.05, Hmin=0.1, NumberOfSamples = 3, qthreshold = 0.1)
# spatial locations: HH, HS, SS. Plot freqcies distribution based on trimming the flanks
OutFLANKResultsPlotter(OF_out, withOutliers = T, NoCorr = T, Hmin = 0.1, binwidth = 0.005, titletext = "Scan for local selection")
# find out what the outliers are
outliers <- which(OF_out$results$OutlierFlag=="TRUE")
outliers

# We can extract info about the outliers by reading in the vcf file and looking at the annotations
vcf1 <- read.vcfR("SSW_all_biallelic.MAF0.02.Miss0.8.recode.vcf")
dim(vcf1)
vcfann <- as.data.frame(getFIX(vcf1))
vcfann[outliers,]
