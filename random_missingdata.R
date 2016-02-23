## Code to insert missing data into our datasets 

setwd("/Users/AmandaXuereb/Documents/PhD/LG_course_project_2016/Simulations/Gsims/")       # change to your directory where the data is stored 

#################
### Data Prep ###
#################

dat <- read.csv("H5S05D15-R5.csv")  # Genetic data; H5 landscape, 5% selection strength, 15% dispersal, one replicate 

# Sample the full data set:
samp <- read.csv("/Users/AmandaXuereb/Documents/PhD/LG_course_project_2016/sample500.csv", header=F)  # 500 sampled individuals
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

#######################################
#######################################

## This is code to loop missing data into all 10 replicates within a directory 


# first, read in the sample500.csv, or sample100.csv, and change to a vector 

samp <- read.csv("/Users/AmandaXuereb/Documents/PhD/LG_course_project_2016/sample500.csv", header=F)  # 500 sampled individuals
samp <- samp[,1]

# set working directory to the directory where replicates are stored; change this to you own path!

setwd("/Users/AmandaXuereb/Documents/PhD/LG_course_project_2016/Simulations/Gsims/GS10D15")

list.dir <- dir() # this will be the 10 replicate files 
list.dir

#read in all the replicate files: 
data <- lapply(list.dir, read.csv)
loci <- seq(from=1, to=200, by=2)

for (i in 1:length(data)) {
  data[[i]] <- data[[i]][samp,]   # sample to 500 individuals (or 100 if samp = sample100.csv)
  data[[i]] <- data[[i]][,8:207]  # remove the columns that are not needed
  data[[i]] <- data[[i]][,loci]   # convert to 1-allele dataset
}

class(data)
dim(data[[1]])
View(data[[1]])

snps <- data  ## just to call it something other than data

# convert each element in the list to a matrix: 
snp.mat <- lapply(snps, function(x) as.matrix(x))  

allsnps_5per <- snp.mat     # make a copy; you can keep the original snp.mat and change the name of the new object to 5per, 10per, 20 per, etc, for percentage of missing data 

# set the amount of missing data to randomly sample:
n <- round(0.05*length(snp.mat[[1]]))   
n

for (j in 1:length(allsnps_5per)) {
  allsnps_5per[[j]][sample(1:length(allsnps_5per[[j]]), n, replace = FALSE)] <- -9
}

head(allsnps_5per[[1]])
lapply(allsnps_5per, function(x) length(which(x == -9)))  # to check that the correct amount of missing data is generated; the length of each list element should be the percentage you are generating * 50,000

names(allsnps_5per) <- seq(0,9,1)
names(allsnps_5per)

sapply(names(allsnps_5per), function(x) write.csv(allsnps_5per[[x]], file = paste(x, "missing5per.csv",sep="_")))







