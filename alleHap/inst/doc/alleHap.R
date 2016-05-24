## ----include=FALSE-------------------------------------------------------
library(knitr)
opts_chunk$set(concordance=TRUE)
library(alleHap)
library(abind)
library(tools)
set.seed(12345678)

## ----echo=FALSE, comment=NA,highlight=TRUE-------------------------------
famIDs <- data.frame(famID="FAM001",indID=1:5,patID=c(0,0,1,1,1),
                     matID=c(0,0,2,2,2),sex=c(1,2,1,2,1),phenot=0)
Mkrs <- rbind(c(1,2, NA,NA, 1,2),c(3,4, 1,2, 3,4),
              c(1,3, 1,2, 1,3),c(NA,NA, 1,1, 2,4),c(1,4, 1,1, 2,4))
colnames(Mkrs)=c("Mk1_1","Mk1_2","Mk2_1","Mk2_2","Mk3_1","Mk3_2")
(family <- cbind(famIDs,Mkrs))

## ----comment=NA,prompt=TRUE,background='#EFF5FB'-------------------------
simulatedFam1 <- alleSimulator(1,2,3,missParProb=0.2,ungenotPars=0.2)
simulatedFam1[[1]]    # Alleles (genotypes) of the 1st simulated family
simulatedFam1[[2]]    # Haplotypes of the 1st simulated family

## ----comment=NA,prompt=TRUE,background='#EFF5FB'-------------------------
simulatedFam2 <- alleSimulator(1,3,3,missOffProb=0.2,ungenotOffs=0.2)
simulatedFam2[[1]]    # Alleles (genotypes) of the 2nd simulated family
simulatedFam2[[2]]    # Haplotypes of the 2nd simulated family

## ----comment=NA,prompt=TRUE,background='#EFF5FB'-------------------------
example1 <- file.path(find.package("alleHap"), "examples", "example1.ped")

example1Alls <- alleLoader(example1)  # Loaded alleles of the example 1
example1Alls[1:10,1:20]  # Alleles of the first 10 subjects

## ----comment=NA,prompt=TRUE,background='#EFF5FB'-------------------------
example2 <- file.path(find.package("alleHap"), "examples", "example2.ped")

example2Alls <- alleLoader(example2)  # Loaded alleles of the example 2
example2Alls[1:9,]  # Alleles of the first 9 subjects

## ----comment=NA,prompt=TRUE,background='#EFF5FB'-------------------------
## Simulation of a family contatining parental missing data
simulatedFam1 <- alleSimulator(1,2,3,missParProb=0.6,dataSummary=FALSE)

simulatedFam1[[1]]  # Alleles of the simulated family

## Allele imputation of the previous family
imputedFam1 <- alleImputer(simulatedFam1[[1]])  
imputedFam1$imputedMkrs # Imputed alleles (markers)

## ----comment=NA,prompt=TRUE,background='#EFF5FB'-------------------------
## Simulation of two families containing offspring missing data
simulatedFam2 <- alleSimulator(2,2,3,missOffProb=0.6,dataSummary=FALSE)
simulatedFam2[[1]] # Alleles of the simulated families

## Allele imputation of the previous familes
imputedFam2 <- alleImputer(simulatedFam2[[1]])
imputedFam2$imputedMkrs # Imputed alleles (markers)

## ----comment=NA,prompt=TRUE,background='#EFF5FB'-------------------------
## Simulation of familes containing parental missing data
simulatedFams1 <- alleSimulator(2,2,6,missParProb=0.2,
                                ungenotPars=0.4,dataSummary=FALSE)
## Reconstruction of previous simulated families
phasedFams1 <- hapPhaser(simulatedFams1[[1]])
simulatedFams1[[1]][,-(1:6)] # Simulated alleles
phasedFams1$phasedMkrs[,-(1:6)] # Imputed/Phased markers
phasedFams1$haplotypes # Phased haplotypes

## ----comment=NA,prompt=TRUE,background='#EFF5FB'-------------------------
## Simulation of families containing offspring missing data
simulatedFams2 <- alleSimulator(2,2,6,missOffProb=0.4,
                                ungenotOffs=0.2,dataSummary=FALSE)
simulatedFams2[[1]][,-(1:6)] # Simulated alleles
## Reconstruction of previous simulated families
phasedFams2 <- hapPhaser(simulatedFams2[[1]],dataSummary=FALSE)
phasedFams2$phasedMkrs[,-(1:6)] # Imputed/Phased markers 
phasedFams2$haplotypes # Phased haplotypes

