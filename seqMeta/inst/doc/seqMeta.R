### R code from vignette source 'seqMeta.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: seqMeta.Rnw:200-214
###################################################
rm(list=ls())
library(seqMeta)
######### load example data:
# contains SNPInfo, phenotyes (pheno1, pheno2)
# genotypes (Z1,Z2), and pedigree information (kins) for pheno2
data(seqMetaExample)
ls()

#Perform study-level analysis:
c1 <- prepScores(Z1, y~bmi+sex, SNPInfo = SNPInfo, data = pheno1)

###save the output, which can be passed to a central location.
study1.out.file <- tempfile()
save(c1, file = study1.out.file)


###################################################
### code chunk number 2: seqMeta.Rnw:226-230
###################################################
c2 <- prepScores(Z2, y~bmi+sex, SNPInfo = SNPInfo, kins = kins, 
                    data = pheno2)
study2.out.file <- tempfile()
save(c2, file = study2.out.file)


###################################################
### code chunk number 3: seqMeta.Rnw:291-297
###################################################
load(study1.out.file)
load(study2.out.file)

meta.results <- skatMeta(c1, c2, SNPInfo = SNPInfo)

head(meta.results)


###################################################
### code chunk number 4: seqMeta.Rnw:321-322
###################################################
study1.results <- skatMeta(c1, SNPInfo = SNPInfo)


###################################################
### code chunk number 5: seqMeta.Rnw:343-345
###################################################
meta.t1.results <- burdenMeta(c1, c2, wts = function(maf){maf < 0.01}, 
                              SNPInfo = SNPInfo)


###################################################
### code chunk number 6: seqMeta.Rnw:348-350
###################################################
meta.t1.results <- burdenMeta(c1, c2, wts =1, 
                              mafRange = c(0,0.01), SNPInfo = SNPInfo)


###################################################
### code chunk number 7: seqMeta.Rnw:353-355
###################################################
meta.mb.results <- burdenMeta(c1, c2, wts = function(maf){1/(maf*(1-maf))}, 
                              SNPInfo = SNPInfo)


###################################################
### code chunk number 8: seqMeta.Rnw:358-359
###################################################
format(head(meta.t1.results),digits=2)


###################################################
### code chunk number 9: seqMeta.Rnw:426-430
###################################################
meta.skato.results <- skatOMeta(c1, c2, rho=seq(0,1,length=11),
     burden.wts = function(maf){dbeta(maf,1,25)}, SNPInfo = SNPInfo, method = "int")

format(head(meta.skato.results),digits=2)


###################################################
### code chunk number 10: seqMeta.Rnw:448-449
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
### code chunk number 13: seqMeta.Rnw:506-509
###################################################
meta.ss.results <- singlesnpMeta(c1, c2, SNPInfo = SNPInfo,
                                 studyBetas = TRUE)
format(head(meta.ss.results),digits=2)


###################################################
### code chunk number 14: seqMeta.Rnw:526-528
###################################################
c1.bin <- prepScores(Z2, ybin~bmi+sex, family = binomial(), 
                     SNPInfo = SNPInfo, data = pheno1)


###################################################
### code chunk number 15: seqMeta.Rnw:534-536
###################################################
c1.cox <- prepCox(Z=Z1, Surv(time, status) ~ bmi + strata(sex), 
                                 SNPInfo = SNPInfo, data =pheno1)


###################################################
### code chunk number 16: seqMeta.Rnw:541-545
###################################################
study1.bin.results <- skatMeta(c1.bin, SNPInfo = SNPInfo, 
                         aggregateBy = "gene")
study1.cox.results <- skatMeta(c1.cox, SNPInfo = SNPInfo, 
                         aggregateBy = "gene")


###################################################
### code chunk number 17: seqMeta.Rnw:562-565
###################################################
meta.results <- skatMeta(c1, c2, SNPInfo = SNPInfo, 
                         aggregateBy = "gene", mafRange = c(0,0.05))
head(meta.results)


###################################################
### code chunk number 18: seqMeta.Rnw:584-599
###################################################
adjustments <- SNPInfo[c(1:3, 20,100), ]
adjustments

####run on each study:
c1.adj <- prepCondScores(Z=Z1, y~sex+bmi, SNPInfo = SNPInfo, 
        adjustments=adjustments, data =pheno1)
c2.adj <- prepCondScores(Z=Z2, y~sex+bmi, SNPInfo = SNPInfo, 
        adjustments=adjustments, kins=kins, data=pheno2)

SNPInfo.sub <- subset(SNPInfo, (SNPInfo$gene %in% adjustments$gene) & 
        !(SNPInfo$Name %in% adjustments$Name) )

#skat
out.skat <- skatMeta(c1.adj,c2.adj, SNPInfo = SNPInfo.sub)
head(out.skat)


###################################################
### code chunk number 19: seqMeta.Rnw:623-638
###################################################
##subset SNPInfo file to first 50 genes, and second 50 genes:
SNPInfo1 <- subset(SNPInfo, SNPInfo$gene %in% unique(SNPInfo$gene)[1:50] )
SNPInfo2 <- subset(SNPInfo, !(SNPInfo$gene %in% unique(SNPInfo$gene)[1:50]) )

##subset corresponding genotype files:
Z1.1 <- subset(Z1, select = colnames(Z1) %in% SNPInfo1$Name) 
Z1.2 <- subset(Z1, select = colnames(Z1) %in% SNPInfo2$Name) 

##run prepScores separately on each chunk:
c1.1 <- prepScores(Z1.1, y~bmi+sex, SNPInfo = SNPInfo1, data = pheno1)
c1.2 <- prepScores(Z1.2, y~bmi+sex, SNPInfo = SNPInfo2, data = pheno1)

##combine results:
c1 <- c(c1.1, c1.2)
class(c1)


