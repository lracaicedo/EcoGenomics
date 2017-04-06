## Laura Caicedo-Quiroga
# Homework 3 for PBIO381

### Analysis of population diversity and structure
# Load the libraries
library(vcfR)
library(adegenet)

#Read the vcf SNP data into R for filter 1
vcf1 <- read.vcfR("SSW_biallelic.MAF0.02.Miss0.7.recode.vcf")

# The adegenet package uses a highly efficient way of storing large SNP datasets in R called a "genlight" object. The following function creates a genlight object from your vcf:
gl1 <- vcfR2genlight(vcf1)
print(gl1) # Looks good! Right # of SNPs and individuals!

# For info, try:
gl1$ind.names #names in the order of vcf files
gl1$loc.names[1:10] #first 10 locus names
gl1$chromosome[1:3] #just transcript IDs

# Notice there's nothing in the field that says "pop"? Let's fix that...
ssw_meta <- read.table("ssw_healthloc.txt", header=T) # read in the metadata
ssw_meta <- ssw_meta[order(ssw_meta$Individual),] # sort by Individual ID, just like the VCF file
str(ssw_meta)

# Confirm the ID's are ordered the same in gl1 and ssw_meta:
gl1$ind.names #see names and sorted from small to big
ssw_meta$Individual

gl1$pop <- ssw_meta$Location # assign locality info
gl1$other <- as.list(ssw_meta$Trajectory)

###DAPC: discriminant analysis on PCA, give a priori groups and ask how well the SNPS differentiate the groups
# Run the DAPC using disease status to group samples
disease.dapc <- dapc(gl1, pop=as.factor(unlist(gl1$other)), n.pca=8, n.da=3,
                     var.loadings=T, pca.info=T)
disease.dapc
# Scatterplot of results
scatter.dapc(disease.dapc, grp=as.factor(unlist(gl1$other)), legend=T)

# Plot the posterior assignment probabilities to each group
compoplot(disease.dapc)
disease.dapc

# Which loci contribute the most to distinguishing Healthy vs. Sick individuals?
loadingplot(abs(disease.dapc$var.load), 
            lab.jitter=1, 
            threshold=quantile(abs(disease.dapc$var.load), probs=0.999))

# Filter #2
#Read the vcf SNP data into R for filter 3
vcf1 <- read.vcfR("SSW_biallelic.MAF0.05.Miss0.7.recode.vcf")

# The adegenet package uses a highly efficient way of storing large SNP datasets in R called a "genlight" object. The following function creates a genlight object from your vcf:
gl1 <- vcfR2genlight(vcf1)
print(gl1) # Looks good! Right # of SNPs and individuals!

# For info, try:
gl1$ind.names #names in the order of vcf files
gl1$loc.names[1:10] #first 10 locus names
gl1$chromosome[1:3] #just transcript IDs

# Notice there's nothing in the field that says "pop"? Let's fix that...
ssw_meta <- read.table("ssw_healthloc.txt", header=T) # read in the metadata
ssw_meta <- ssw_meta[order(ssw_meta$Individual),] # sort by Individual ID, just like the VCF file
str(ssw_meta)

# Confirm the ID's are ordered the same in gl1 and ssw_meta:
gl1$ind.names #see names and sorted from small to big
ssw_meta$Individual

gl1$pop <- ssw_meta$Location # assign locality info
gl1$other <- as.list(ssw_meta$Trajectory)

###DAPC: discriminant analysis on PCA, give a priori groups and ask how well the SNPS differentiate the groups
# Run the DAPC using disease status to group samples
disease.dapc <- dapc(gl1, pop=as.factor(unlist(gl1$other)), n.pca=8, n.da=3,
                     var.loadings=T, pca.info=T)
disease.dapc
# Scatterplot of results
scatter.dapc(disease.dapc, grp=as.factor(unlist(gl1$other)), legend=T)

# Plot the posterior assignment probabilities to each group
compoplot(disease.dapc)
disease.dapc

# Which loci contribute the most to distinguishing Healthy vs. Sick individuals?
loadingplot(abs(disease.dapc$var.load), 
            lab.jitter=1, 
            threshold=quantile(abs(disease.dapc$var.load), probs=0.999))

