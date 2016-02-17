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

############################
## Missing adaptive locus ##
############################

## This is code for creating missing data for the adaptive locus in a certain percentage of the individuals 
## Missing data is represented as -9 

## Test small fake dataset: 

#ind_a <- c(0,0,1,2,1,1)
#ind_b <- c(1,1,1,2,2,0)
#ind_c <- c(2,1,0,1,0,0)
#ind_d <- c(2,2,1,1,0,1)
#ind_e <- c(0,1,1,2,2,0)
#ind_f <- c(1,1,1,1,0,1)

#test_snps<-rbind(ind_a,ind_b,ind_c,ind_d,ind_e,ind_f)
#colnames(test_snps) <- paste0("loc_", seq(1,6,1))
#test_snps # loc_6 will be the adaptive locus 

# replace 50% of locus 6 with -9

#n <- round(0.5*length(test_snps[,6]))
#n

#test_snps_loc6 <- test_snps[,6]
#test_snps_loc6

#test_snps_loc6[sample(1:length(test_snps_loc6), n, replace=FALSE)] <- -9
#test_snps_loc6

#test_snps_50per <- test_snps
#test_snps_50per[,6] <- test_snps_loc6
#test_snps_50per


### Real SNP data ###

dim(snps)  # 100 x 500

n <- round(0.1*length(snps[,100])) # 10% of individuals (rows); change for other percentages 
n

sel_snp <- snps[,1]
sel_snp
sel_snp == snps[,1] # to check that they are equal to each other; this should all be TRUE 

sel_snp[sample(1:length(sel_snp), n, replace=FALSE)] <- -9
sel_snp
length(which(sel_snp == -9)) # check how many -9s there are 


snps_sel10per <- snps       # make a copy of the snps dataset
snps_sel10per[,1] <- sel_snp
head(snps_sel10per)
length(which(snps_sel10per[,1] == -9)) # check how many -9s there are in column 100
which(snps_sel10per[,2:100] == -9) # check to make sure there are no -9s in the other columns 

##########################
## Missing neutral loci ##
##########################

head(snps)        
length(snps)

snps_neutral <- snps[,-1] # remove the adaptive locus 
ncol(snps_neutral) 

snps_neut.mat <- as.matrix(snps_neutral)
length(snps_neut.mat)


n <- round(0.5*length(snps_neut.mat))   # change 0.5 (50%) to other percentages of missing data
n 

snps_neut50per <- snps_neut.mat            
snps_neut50per[sample(1:length(snps_neut50per), n, replace=FALSE)] <- -9
head(snps_neut50per) 
snps_neut50per <- as.data.frame(snps_neut50per)

length(which(snps_neut50per == -9))  # to see if the correct number of snps have been changed to -9

###############################################
## Combine missing adaptive and neutral loci ##
###############################################

head(snps_sel10per)
dim(snps_sel10per[,2:100]) == dim(snps_neut50per)

snps_combined_missing <- snps_sel10per # make a copy of the dataset with missing data at the adaptive locus
snps_combined_missing[,2:100] <- snps_neut50per

length(which(snps_combined_missing[,2:100]==-9))
length(which(snps_combined_missing == -9))


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


