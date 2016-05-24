### R code from vignette source 'skatMeta.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: skatMeta.Rnw:202-216
###################################################
rm(list=ls())
library(skatMeta)
######### load example data:
# contains SNPInfo, phenotyes (pheno1, pheno2)
# genotypes (Z1,Z2), and pedigree information (kins) for pheno2
data(skatExample)
ls()

#Perform cohort-level analysis:
c1 <- skatCohort(Z1, y~bmi+sex, SNPInfo = SNPInfo, data = pheno1)

###save the output, which can be passed to a central location.
cohort1.out.file <- tempfile()
save(c1, file = cohort1.out.file)


###################################################
### code chunk number 2: skatMeta.Rnw:233-237
###################################################
c2 <- skatFamCohort(Z2, y~bmi+sex, SNPInfo = SNPInfo, fullkins = kins, 
                    data = pheno2)
cohort2.out.file <- tempfile()
save(c2, file = cohort2.out.file)


###################################################
### code chunk number 3: skatMeta.Rnw:298-304
###################################################
load(cohort1.out.file)
load(cohort2.out.file)

meta.results <- skatMeta(c1, c2, SNPInfo = SNPInfo)

head(meta.results)


###################################################
### code chunk number 4: skatMeta.Rnw:328-329
###################################################
cohort1.results <- skatMeta(c1, SNPInfo = SNPInfo)


###################################################
### code chunk number 5: skatMeta.Rnw:350-352
###################################################
meta.t1.results <- burdenMeta(c1, c2, wts = function(maf){maf < 0.01}, 
                              SNPInfo = SNPInfo)


###################################################
### code chunk number 6: skatMeta.Rnw:355-357
###################################################
meta.t1.results <- burdenMeta(c1, c2, wts =1, 
                              mafRange = c(0,0.01), SNPInfo = SNPInfo)


###################################################
### code chunk number 7: skatMeta.Rnw:360-362
###################################################
meta.mb.results <- burdenMeta(c1, c2, wts = function(maf){1/(maf*(1-maf))}, 
                              SNPInfo = SNPInfo)


###################################################
### code chunk number 8: skatMeta.Rnw:365-366
###################################################
format(head(meta.t1.results),digits=2)


###################################################
### code chunk number 9: skatMeta.Rnw:433-437
###################################################
meta.skato.results <- skatOMeta(c1, c2, rho=seq(0,1,length=11),
     burden.wts = function(maf){dbeta(maf,1,25)}, SNPInfo = SNPInfo, method = "int")

format(head(meta.skato.results),digits=2)


###################################################
### code chunk number 10: skatMeta.Rnw:455-456
###################################################
table(meta.skato.results$rho)


###################################################
### code chunk number 11: fig1plot
###################################################
wu.burden <- burdenMeta(c1, c2, wts = function(maf){dbeta(maf,1,25)}, 
                        SNPInfo=SNPInfo)
pseq <- seq(0,1,length=100)
plot(y=meta.skato.results$p, x=pmin(wu.burden$p,meta.results$p), 
     xlab ="Minimum of SKAT and Burden", ylab = "SKAT-O")
abline(0,1,lty=2)
lines(x=pseq,y=1-(1-pseq)^2,col=2,lty=2,lwd=2)
legend("bottomright", lwd=2,lty=2,col=2,legend="Sidak correction")	


###################################################
### code chunk number 12: fig1
###################################################
wu.burden <- burdenMeta(c1, c2, wts = function(maf){dbeta(maf,1,25)}, 
                        SNPInfo=SNPInfo)
pseq <- seq(0,1,length=100)
plot(y=meta.skato.results$p, x=pmin(wu.burden$p,meta.results$p), 
     xlab ="Minimum of SKAT and Burden", ylab = "SKAT-O")
abline(0,1,lty=2)
lines(x=pseq,y=1-(1-pseq)^2,col=2,lty=2,lwd=2)
legend("bottomright", lwd=2,lty=2,col=2,legend="Sidak correction")	


###################################################
### code chunk number 13: skatMeta.Rnw:513-516
###################################################
meta.ss.results <- singlesnpMeta(c1, c2, SNPInfo = SNPInfo,
                                 cohortBetas = TRUE)
format(head(meta.ss.results),digits=2)


###################################################
### code chunk number 14: skatMeta.Rnw:533-535
###################################################
c1.bin <- skatCohort(Z2, ybin~bmi+sex, family = binomial(), 
                     SNPInfo = SNPInfo, data = pheno1)


###################################################
### code chunk number 15: skatMeta.Rnw:541-543
###################################################
c1.cox <- skatCoxCohort(Z=Z1, Surv(time, status) ~ bmi + strata(sex), 
                                 SNPInfo = SNPInfo, data =pheno1)


###################################################
### code chunk number 16: skatMeta.Rnw:548-552
###################################################
cohort1.bin.results <- skatMeta(c1.bin, SNPInfo = SNPInfo, 
                         aggregateBy = "gene")
cohort1.cox.results <- skatMeta(c1.cox, SNPInfo = SNPInfo, 
                         aggregateBy = "gene")


###################################################
### code chunk number 17: skatMeta.Rnw:569-572
###################################################
meta.results <- skatMeta(c1, c2, SNPInfo = SNPInfo, 
                         aggregateBy = "gene", mafRange = c(0,0.05))
head(meta.results)


###################################################
### code chunk number 18: skatMeta.Rnw:591-606
###################################################
adjustments <- SNPInfo[c(1:3, 20,100), ]
adjustments

####run on each cohort:
c1.adj <- skatCohortAdj(Z=Z1, y~sex+bmi, SNPInfo = SNPInfo, 
        adjustments=adjustments, data =pheno1)
c2.adj <- skatFamCohortAdj(Z=Z2, y~sex+bmi, SNPInfo = SNPInfo, 
        adjustments=adjustments, fullkins=kins, data=pheno2)

SNPInfo.sub <- subset(SNPInfo, (SNPInfo$gene %in% adjustments$gene) & 
        !(SNPInfo$Name %in% adjustments$Name) )

#skat
out.skat <- skatMeta(c1.adj,c2.adj, SNPInfo = SNPInfo.sub)
head(out.skat)


###################################################
### code chunk number 19: skatMeta.Rnw:630-645
###################################################
##subset SNPInfo file to first 50 genes, and second 50 genes:
SNPInfo1 <- subset(SNPInfo, SNPInfo$gene %in% unique(SNPInfo$gene)[1:50] )
SNPInfo2 <- subset(SNPInfo, !(SNPInfo$gene %in% unique(SNPInfo$gene)[1:50]) )

##subset corresponding genotype files:
Z1.1 <- subset(Z1, select = colnames(Z1) %in% SNPInfo1$Name) 
Z1.2 <- subset(Z1, select = colnames(Z1) %in% SNPInfo2$Name) 

##run skatCohort separately on each chunk:
c1.1 <- skatCohort(Z1.1, y~bmi+sex, SNPInfo = SNPInfo1, data = pheno1)
c1.2 <- skatCohort(Z1.2, y~bmi+sex, SNPInfo = SNPInfo2, data = pheno1)

##combine results:
c1 <- c(c1.1, c1.2)
class(c1)


