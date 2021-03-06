## This code will loop through all replicates for a set of 24 simulations (Gsims, H1sims, H5sims, H9sims) and 6 simulations for Nsims, with 1 missing data input file.  Keep the same file structure as when you download the files from Dryad (i.e. one folder per sim, and 10 replicate files in each folder named "R0.csv" to "R9.csv"). 


###########
## TO DO ##
###########

##### H9 SIMS
  ###  500 individuals  
   ##   neutral missing:  5% 10% 20% 50%
   ##   adaptive missing: 5% 10% 20% 50%
   ##   no missing 
  ###  100 individuals 
   ##   random missing:   5% 10% 20% 50%
   ##   neutral missing:  5% 10% 20% 50%
   ##   adaptive missing: 5% 10% 20% 50%
   ##   no missing   

##### H5 SIMS
  ###  500 individuals 
   ##   random missing:   5% 10% 20% 50%
   ##   neutral missing:  5% 10% 20% 50%
   ##   adaptive missing: 5% 10% 20% 50%
   ##   no missing 

   ### 100 individuals 
    ##  random missing:   5% 10% 20% 50%
    ##  neutral missing:  5% 10% 20% 50%
    ##  adaptive missing: 5% 10% 20% 50%
    ##  no missing   

##### H1 SIMS
  ###  500 individuals 
   ##   random missing:   5% 10% 20% 50%
   ##   neutral missing:  5% 10% 20% 50%
   ##   adaptive missing: 5% 10% 20% 50%
   ##   no missing 

  ### 100 individuals 
   ##  random missing:   5% 10% 20% 50%
   ##  neutral missing:  5% 10% 20% 50%
   ##  adaptive missing: 5% 10% 20% 50%
   ##  no missing  

##### G SIMS
  ###  500 individuals 
   ##   random missing:   5% 10% 20% 50%
   ##   neutral missing:  5% 10% 20% 50%
   ##   adaptive missing: 5% 10% 20% 50%
   ##   no missing 

  ### 100 individuals 
   ##  random missing:   5% 10% 20% 50%
   ##  neutral missing:  5% 10% 20% 50%
   ##  adaptive missing: 5% 10% 20% 50%
   ##  no missing  

##########
## DONE ##
##########

##### H9 SIMS
###  500 individuals  
##   random missing:   5% 10% 20% 50%

rm(list=ls())
library(raster)
library(vegan)

###set working directory to the folder that stores the 240 sim folders that you are running
setwd("/Users/AmandaXuereb/Documents/PhD/LG_course_project_2016/Simulations/H9sims")
sim <- dir()#[-c(1,15,18,21)]#for H5 sims 

repl <- c(0:9)


##########################
### Loop starts here 
##########################

for (s in 1:length(sim)) {
  for (r in 1:length(repl)) {
    
###read in the habitat files
###R1 goes with sim rep0
###Change H9, H5, H1
    fname <- paste("/Users/AmandaXuereb/Documents/PhD/LG_course_project_2016/Surfaces_Sample/L10H9R",r,"_aa.asc", sep="")
    qrule <- raster(fname)                                

    
###read in the snp dataset
###Change directory for H9, H5, H1
    fname1 <- paste("/Users/AmandaXuereb/Documents/PhD/LG_course_project_2016/Simulations/H9sims/",sim[s],"/R",repl[r],".csv", sep="")
    data <- read.csv(fname1)

###extract coordinates and the habitat value at each coordinate    
    coord <- data[,2:3]
    habitat <- extract(qrule,coord)

###convert UTM - subtract corner coordinate and /1000 to convert from m to km
    xvar <- (data[,2] - 477895.2)/1000         
    yvar <- (data[,3] - 4905702)/1000 

###create the full data set
    fulldata <- cbind(xvar, yvar, habitat, data[,8:207])
    
###subset individuals from the sample100 or sample500.csv files
    samp <- read.csv("/Users/AmandaXuereb/Documents/PhD/LG_course_project_2016/sample30.csv", header=F)
    samp <- samp[,1]  
    sampdata <- fulldata[samp,]

###extract the environmental variables alone  
    env <- sampdata[,1:3] 
    env <- scale(env, center=T, scale=T)
    env <- data.matrix(env)
    colnames(env)<-c("xvar","yvar","habitat")
    
###extract the snps alone
    snps <- sampdata[,4:203]
    snps <- snps[,seq(1, ncol(snps), 2)]
    colnames(snps) <- seq(0, 99, by=1)
    snps <- data.matrix(snps)
    
#enter NAs
##in glob2rx(), change "random_", "neutral_", or "adaptive" and "_100" or "_500"   
    missing_files <- list.files("/Users/AmandaXuereb/Documents/PhD/LG_course_project_2016/Simulations/Input Files/", pattern = glob2rx("neutral_*_30.csv"), full.names = F, recursive = F) 
    
###change missing_files[x] to [1] to [4] for each percentage of missing data
    missing <- read.csv(paste0("/Users/AmandaXuereb/Documents/PhD/LG_course_project_2016/Simulations/Input Files/", missing_files[1]),header=FALSE)  
    missing <- as.numeric(missing[,1])  
    snps[missing] <- NA
    
    
###remove NA rows (individuals) from both snp and env datasets; there are a few replicates with NA individuals so there will be 499 instead of 500 individuals
    subsetNA <- apply(snps,1,function(x) length(unique(x))>1)
    snps <- snps[subsetNA,]                                      
    env <- env[subsetNA,] 
    

###remove monomorphic snps (I don't think we have any, but just in case)
    subsetMM <- apply(snps,2,function(x) length(unique(x))>1)
    snps <- snps[,subsetMM]                                     
    
###MAF filtering 
    geno <- as.matrix(snps)
    
## calc_n
    n0 <- apply(geno==0,2,sum,na.rm=T)
    n1 <- apply(geno==1,2,sum,na.rm=T)
    n2 <- apply(geno==2,2,sum,na.rm=T)
    n <- n0 + n1 + n2
    
## calculate allele frequencies
    p <- ((2*n0)+n1)/(2*n)
    q <- 1 - p
    maf <- pmin(p, q)
    
## KEEP LOCI WITH MAF > 3%
    subsetMAF <- maf[maf>0.03]
    tempname <- maf[maf<0.03]
    nameMAF <- match(names(subsetMAF), colnames(snps))
    newsnps <- snps
    snps <- snps[, c(nameMAF)]
    
    
    
###To check which loci have been removed & final # individuals, this will create output files further down so we can just see how many individuals/snps are remaining after removing NAinds, maf < 3%, and to double check the number of NA snps
    #tempname <- names(subsetMM)
    dims <- dim(snps)
    ind_missing <- cbind(sim[s], repl[r], dims[1], dims[2], length(tempname))
    fname <- paste("NAinds_", sim[s], repl[r], sep="")
    assign(fname, ind_missing)
    
#     tempname <- names(subsetMAF)
#     dims <- dim(snps)
#     maf <- cbind(sim[s], repl[r], dims[1], dims[2], tempname[subsetMAF==F])
#     fname <- paste("maf_", sim[s], repl[r], sep="")
#     assign(fname, maf)
    
    snps_missing <- cbind(sim[s],repl[r],length(which(is.na(snps)==TRUE)))
    fname <- paste("NAsnps_", sim[s], repl[r], sep="")
    assign(fname, snps_missing)
    

## You shouldn't have to change anything inside the loop from here; but you will have to change the output file names outside of the loop 

    ##Impute missing data with colMeans
    
    #colMeans(snps, na.rm = TRUE)
    #mean(snps[,1], na.rm=TRUE)
    
    cM <- colMeans(snps, na.rm =TRUE)
    indx <- which(is.na(snps), arr.ind=TRUE)
    snps[indx] <- cM[indx[,2]]
    
    cM2 <- colMeans(newsnps, na.rm=TRUE)
    indx2 <- which(is.na(newsnps), arr.ind=TRUE)
    newsnps[indx2] <- cM2[indx2[,2]]
   
    
###format the snp data
    snps.scale <- scale(snps, center=T, scale=T)  ## scale and center for PCA/RDA
    snps.bray <- vegdist(snps, method="bray")  ## bray-curtis distance for PCoA and dbRDA
    
############
##  RDA  ##
#############
    snp.rda <- rda(snps.scale, env, scale=F)
    snp.rda.sum <- summary(snp.rda)
    snpload.rda <- snp.rda.sum$species[,1:3]
    
##############
##  dbRDA  ##
##############
    snp.dbrda <- capscale(snps.bray ~ env, comm=snps)  # RDA on a PCoA
    snp.dbrda.sum <- summary(snp.dbrda)
    snpload.dbrda <- snp.dbrda.sum$species[,1:3]   # just like in RDA
    rownames(snpload.dbrda) <- colnames(snps)
    
#####################    
## detect outliers ##
#####################    
    tests <- list(snpload.rda, snpload.dbrda)
    
    test.names <- c("rda","dbrda")
    
    for (i in 1:length(tests)) {
      
      for (j in 1:3) {                                    
        x <- tests[i] #change back to i 
        x <- matrix(unlist(x), ncol = 3)
        rownames(x) <- colnames(snps)
        x <- x[,j] # change back to j 
        lims <- mean(x) + c(-1, 1) * 3 * sd(x)         
        out <- x[x < lims[1] | x > lims[2]]
        
        if (length(out) > 0) {
          names <- as.numeric(names(out))                    
          axis <- rep(j, length(out))
          outdata <- t(rbind(axis, names, as.numeric(out)))
          
          snpcors <- matrix(NA, nrow=length(names), ncol=6)
          
          for (k in 1:length(names)) {
            outlocus <- names[k] + 1   # to access correct column in snp dataframe
            
            if (ncol(snps) == 100) { 
              outsnp <- snps[,outlocus]
            } else if (ncol(snps) < 100) {
              outsnp <- newsnps[,outlocus] # use the full snp matrix including the MAF < 3% locus 
              outsnp <- outsnp[!is.na(outsnp)] # remove NA individuals, if any
            }  
            
            corr.x <- cor.test(outsnp, env[,1])
            snpcors[k,1] <- corr.x$estimate
            snpcors[k,2] <- corr.x$p.value
            
            corr.y <- cor.test(outsnp, env[,2])
            snpcors[k,3] <- corr.y$estimate
            snpcors[k,4] <- corr.y$p.value
            
            corr.h <- cor.test(outsnp, env[,3])
            snpcors[k,5] <- corr.h$estimate
            snpcors[k,6] <- corr.h$p.value
          }
          
          outdata <- cbind(outdata, snpcors)
          fname <- paste(test.names[i],j, sep="")
          assign(fname, outdata)
        }
        else if (length(out) == 0) {
          fname <- paste(test.names[i],j, sep="")
          assign(fname, NA)
        }
      }    
      
    }
    
    out.rda <- rbind(rda1, rda2, rda3)
    out.rda <- as.data.frame(out.rda)
    label0 <- rep("RDA", nrow(out.rda))
    label1 <- rep(sim[s], nrow(out.rda))
    label2 <- rep(repl[r], nrow(out.rda))
    out.rda <- cbind(label0,label1,label2,out.rda)
    out.rda <- out.rda[complete.cases(out.rda),]
    
    out.dbrda <- rbind(dbrda1, dbrda2, dbrda3)
    out.dbrda <- as.data.frame(out.dbrda)
    label0 <- rep("dbRDA", nrow(out.dbrda))
    label1 <- rep(sim[s], nrow(out.dbrda)) 
    label2 <- rep(repl[r], nrow(out.dbrda))
    out.dbrda <- cbind(label0,label1,label2,out.dbrda)
    out.dbrda <- out.dbrda[complete.cases(out.dbrda),]
    
    outs <- rbind(out.rda, out.dbrda)
    
    if (nrow(outs) > 0) {   
      colnames(outs) <- c("ord","sim","rep", "axis","locus","loading", "x-corr", "x-pval", "y-corr", "y-pval", "h-corr", "h-pval")
      fname <- paste("out_", sim[s],"_R",repl[r], sep="")
      assign(fname,outs)  
    }
    
  }
}

########################
##### Loop ends here
##########################


## save output; there will be 5 output files in total
# 1. Raw Data
# 2. True Positives 
# 3. Summary file 
# 4. Number of individuals after removing NA individuals 
# 5. Number of snps converted to NAs 

## Follow this output file naming : "_neutral_50percent_500ind_RawData.csv" with the conditions that you ran

save <- ls()[grep("out_", ls())]

bar = NA

for (l in 1:length(save)) {
  foo <- get(save[l])
  bar <- rbind(foo, bar)
}

bar <- bar[-nrow(bar),]

### CHANGE OUTPUT FILE NAME: 
fname <- paste("/Users/AmandaXuereb/Documents/PhD/LG_course_project_2016/Simulations/H9sims_output/H9_neutral/", sim[s], "_neutral_5percent_30ind_imputelast_RawData.csv", sep="")
write.csv(bar, file=fname)


## true positives:
tp <- bar[bar$locus == 0,]

### CHANGE OUTPUT FILE NAME: 
fname <- paste("/Users/AmandaXuereb/Documents/PhD/LG_course_project_2016/Simulations/H9sims_output/H9_neutral/", sim[s], "_neutral_0percent_30ind_TruePos.csv", sep="")
write.csv(tp, file=fname)

# summary of true and false positives
fp <- bar[bar$locus != 0,]
ord <- c("RDA","dbRDA")
summary <- matrix(NA, nrow=length(sim)*2, ncol=5)
colnames(summary) <- c("sim", "ord", "tp", "fp", "fp.sd")

summary[,1] <- as.vector(sapply(sim, function(x) rep(x,2)))
summary[,2] <- as.vector(rep(ord,length(sim)))


for (i in 1:length(sim)) {
  foo <- tp[tp$sim == sim[i],]
  baz <- fp[fp$sim == sim[i],]
  
  for (j in 1:length(ord)) {
    
    # true positives    
    bar <- foo[foo$ord == ord[j],]
    bar1 <- bar[!duplicated(bar$rep),]
    rowindex <- (i-1)*2
    summary[j+rowindex,3] <- nrow(bar1)
    
    #false positives
    qux <- baz[baz$ord == ord[j],]
    
    temp <- vector(mode="integer", length=10)
    
    for (k in 1:length(repl)) {
      rux <- qux[qux$rep == repl[k],]
      rux1 <- rux[!duplicated(rux$locus),] # removes any loci detected on >1 axis in same replicate
      temp[k] <- nrow(rux1)
    }
    
    summary[j+rowindex,4] <- sum(temp)
    summary[j+rowindex,5] <- sd(temp)
  }
}

### CHANGE OUTPUT FILE NAME: 
fname <- paste("/Users/AmandaXuereb/Documents/PhD/LG_course_project_2016/Simulations/H9sims_output/H9_neutral/", sim[s], "_neutral_0percent_30ind_Summary.csv", sep="")
write.csv(summary, file=fname)


## save file telling how many NA individuals removed 

library(gtools)

save <- ls()[grep("NAinds_", ls())]

bar = data.frame()

for (l in 1:length(save)) {
  foo <- get(save[l])
  foo <- as.data.frame(foo)
  bar <- smartbind(foo, bar)
}

bar <- bar[-nrow(bar),]
colnames(bar) <- c("sim","repl","num inds", "loci","num loci MAF < 3%")
rownames(bar) <- seq(1,nrow(bar),1)

### CHANGE OUTPUT FILE NAME: 
fname <- paste("/Users/AmandaXuereb/Documents/PhD/LG_course_project_2016/Simulations/H9sims_output/H9_neutral/", sim[s], "_neutral_0percent_30ind_NAinds.csv", sep="")
write.csv(bar, file=fname)

## save file telling how many snps are NA

save <- ls()[grep("NAsnps_", ls())]

bar = data.frame()

for (l in 1:length(save)) {
  foo <- get(save[l])
  foo <- as.data.frame(foo)
  bar <- smartbind(foo, bar)
}

bar <- bar[-nrow(bar),]
colnames(bar) <- c("sim","repl","num NA snps")
rownames(bar) <- seq(1, nrow(bar), 1)

### CHANGE OUTPUT FILE NAME: 
fname <- paste("/Users/AmandaXuereb/Documents/PhD/LG_course_project_2016/Simulations/H9sims_output/H9_neutral/", sim[s], "_neutral_50percent_30ind_NAsnps.csv", sep="")
write.csv(bar, file=fname)


# save <- ls()[grep("maf_", ls())]
# 
# bar = data.frame()
# 
# for (l in 1:length(save)) {
#   foo <- get(save[l])
#   foo <- as.data.frame(foo)
#   bar <- smartbind(foo, bar)
# }
# 
# bar <- bar[-nrow(bar),]
# colnames(bar) <- c("sim","repl","num NA snps")
# rownames(bar) <- seq(1, nrow(bar), 1)
# 
# ### CHANGE OUTPUT FILE NAME: 
# fname <- paste("/Users/AmandaXuereb/Documents/PhD/LG_course_project_2016/Simulations/H5sims_output/H5_random/", sim[s], "_random_10percent_500ind_maf.csv", sep="")
# write.csv(bar, file=fname)
# 
# 
