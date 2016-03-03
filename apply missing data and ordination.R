
rm(list=ls())

library(raster)
library(vegan)

### I created a folder with all of the datasets and replicates so that the file name of the data set reflects the habitat, selection strength, dispersal level, and replicate numbers 0-9
### For example H1S01D03_R0.csv ... and they are in a folder called "ALLsims" 

setwd("/Users/AmandaXuereb/Documents/PhD/LG_course_project_2016/Simulations/ALLsims/")
dir()

# create a vector of dataset names 
sims <- list.files(pattern = glob2rx("H1*"),full.names=F, recursive=F)    # creates a list of file names that start with H1 (for example, to start with one habitat type). This will have to be changed for G, H5, H9, and N

sims
length(sims) # 240 (24 * 10 reps)


# create a vector of replicate numbers 
repl <- c(0:9)

# create a vector of missing data scenario input files. I would just stick to using 4 at a time for now (for example random, 4 levels, and n=100); there are 20 input files: 3 scenarios (random, adaptive locus only, neutral loci only) x 4 levels (5, 10, 20, 50 %) x 2 sample sizes (100, 500)

missing <- list.files("/Users/AmandaXuereb/Documents/PhD/LG_course_project_2016/Simulations/", pattern = glob2rx("random_*_100.csv"), full.names = F, recursive = F)  

missing


# this is to select every other column in the snp dataset later on
loci <- seq(from=1, to=200, by=2)      



## Start looping

## you will need:
## sample100 or sample500 .csv file
## "Surfaces_Sample" folder from Dryad 
## missing data .csv input files 


# Start by working with one sample size and one habitat configuration, and re-run for others. 


for (s in 1:length(sims)) {
  for (r in 1:length(repl)) {  # 
    
    fname <- paste("/Users/AmandaXuereb/Documents/PhD/LG_course_project_2016/Surfaces_Sample/L10H1R",r,"_aa.asc", sep="")
    qrule <- raster(fname)                                    # CHANGE FOR H1/H5/H9 -- R1 goes with sim rep0
    
    fname1 <- paste("/Users/AmandaXuereb/Documents/PhD/LG_course_project_2016/Simulations/ALLsims/",sims[s], sep="")
    data <- read.csv(fname1)
    coord <- data[,2:3]
    habitat <- extract(qrule,coord)
    
    xvar <- (data[,2] - 477895.2)/1000         # convert UTM - subtract corner coordinate and /1000 to convert from m to km
    yvar <- (data[,3] - 4905702)/1000 
    
    fulldata <- cbind(xvar, yvar, habitat, data[,8:207])
    
    samp <- read.csv("/Users/AmandaXuereb/Documents/PhD/LG_course_project_2016/sample100.csv", header=F)
    samp <- samp[,1]  
    sampdata <- fulldata[samp,]
    
    env <- sampdata[,1:3] 
    env <- scale(env, center=T, scale=T) # scale and center env vars so they are comparable
    env <- data.matrix(env)
    colnames(env)<-c("xvar","yvar","habitat")
    
    ## fix the snp dataset
    snps <- sampdata[,4:203]
    snps <- snps[,seq(1, ncol(snps), 2)]
    colnames(snps) <- seq(0, 99, by=1)
    snps.mat <- as.matrix(snps)
    
    for (m in 1:length(missing)) {
      fname2 <- paste("/Users/AmandaXuereb/Documents/PhD/LG_course_project_2016/Simulations/",missing[m], sep="")
      missing_snps <- read.csv(fname2, header = FALSE)
      missing_snps <- missing_snps[,1]
      snps.mat[missing_snps] <- 'NA'
      snps.mat[snps.mat != 'NA'] <- as.integer(snps.mat != 'NA')
      snps <- as.data.frame(snps.mat)
    

## STOP HERE FOR NOW -- need to make sure it's looping correctly and think about how to integrate ordinations.

## WE need to impute the missing data for RDA/dbRDA to work
      
## format the snp data
      snps.scale <- scale(snps, center=T, scale=T, na.rm=T)
## scale and center for PCA/RDA
      snps.bray <- vegdist(snps, method="bray")  ## bray-curtis distance for PCoA and dbRDA
      
      ##  RDA  ##
      snp.rda <- rda(snps.scale, env, scale=F)
      snp.rda.sum <- summary(snp.rda)
      snpload.rda <- snp.rda.sum$species[,1:3]
      
      ##  dbRDA  ##
      snp.dbrda <- capscale(snps.bray ~ env, comm=snps)  # RDA on a PCoA
      snp.dbrda.sum <- summary(snp.dbrda)
      snpload.dbrda <- snp.dbrda.sum$species[,1:3]   # just like in RDA
      rownames(snpload.dbrda) <- colnames(snps)
      
    }
  }
 