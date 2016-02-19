## Code to insert missing data into our datasets 

setwd("/path/to/file/")       # change to your directory where the data is stored 

#################
### Data Prep ###
#################

dat <- read.csv("H5S05D15-R5.csv")  # Genetic data; H5 landscape, 5% selection strength, 15% dispersal, one replicate 

# Sample the full data set:
samp <- read.csv("sample500.csv", header=F)  # 500 sampled individuals
samp <- samp[,1]                             # convert to a vector
sampdat <- dat[samp,]                        # subsample the full dataset

# Reduce the SNP data set to only one column per locus (from two)
snps <- sampdat[,8:207] 
loci <- seq(from=1, to=200, by=2)
newsnps <- snps[,loci]
newsnps

snps <- newsnps; colnames(snps) <- seq(0, 99, by=1)

#########################################
## Missing data random across all loci ##
#########################################

head(snps)        # use the new data with one allele per locus (matrix of 500 ind x 100 loci)
length(snps)

snp.mat <- as.matrix(snps)
head(snp.mat)
dim(snp.mat)

n <- round(0.1*length(snp.mat))   # change 0.1 (10%) to other percentages of missing data
n 

allsnps_10per <- snp.mat            # make a copy of the original snp data
allsnps_10per[sample(1:length(allsnps_10per), n, replace=FALSE)] <- -9
head(allsnps_10per) 

length(which(allsnps_10per == -9))  #for 10% missing data, the length should = 5000