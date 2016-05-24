### R code from vignette source 'giant_package_vignette.Snw'
### Encoding: UTF-8

###################################################
### code chunk number 1: giant_package_vignette.Snw:45-46 (eval = FALSE)
###################################################
## install.packages("GiANT")


###################################################
### code chunk number 2: giant_package_vignette.Snw:51-56 (eval = FALSE)
###################################################
## # CRAN packages
## install.packages(c("st", "fdrtool"))
## # Bioconductor packages
## source("http://bioconductor.org/biocLite.R")
## biocLite(c("GlobalAncova", "limma", "DESeq"))


###################################################
### code chunk number 3: giant_package_vignette.Snw:61-62
###################################################
library(GiANT)


###################################################
### code chunk number 4: giant_package_vignette.Snw:80-81
###################################################
set.seed(42)


###################################################
### code chunk number 5: giant_package_vignette.Snw:84-99
###################################################
# load data
require(GlobalAncova)
data(vantVeer)
data(phenodata)
data(pathways)

resGsea <- geneSetAnalysis(
	labs = phenodata$metastases,
	method = "pearson",
	numSamples = 1000,
	dat = vantVeer,
	geneSets = pathways,
	analysis = analysis.gsea(),
	adjustmentMethod = "fdr",
	signLevel=0.1)


###################################################
### code chunk number 6: giant_package_vignette.Snw:104-105
###################################################
summary(resGsea)


###################################################
### code chunk number 7: giant_package_vignette.Snw:109-111
###################################################
tab <- createSummaryTable(resGsea)
tab


###################################################
### code chunk number 8: giant_package_vignette.Snw:114-116
###################################################
tab <- createSummaryTable(resGsea, significantOnly=TRUE, orderBy="geneSetName")
tab


###################################################
### code chunk number 9: giant_package_vignette.Snw:129-130 (eval = FALSE)
###################################################
## hist(resGsea, subset = 3, aggregate = TRUE)


###################################################
### code chunk number 10: giant_package_vignette.Snw:132-135
###################################################
pdf("gsea.pdf")
hist(resGsea, subset = 3, aggregate = TRUE)
dev.off()


###################################################
### code chunk number 11: giant_package_vignette.Snw:142-158 (eval = FALSE)
###################################################
## library(parallel)
## mc <- 2 #number of cpus to use
## cl <- makeCluster(mc) #initialize a cluster
## 
## resGsea <- geneSetAnalysis(
## 	labs = phenodata$metastases,
## 	method = "pearson",
## 	numSamples = 1000,
## 	dat = vantVeer,
## 	geneSets = pathways,
## 	analysis = analysis.gsea(),
## 	adjustmentMethod = "fdr",
## 	signLevel=0.1,
## 	cluster = cl)
## 
## stopCluster(cl)


###################################################
### code chunk number 12: giant_package_vignette.Snw:165-166
###################################################
set.seed(132)


###################################################
### code chunk number 13: giant_package_vignette.Snw:169-181
###################################################
stat <- abs(apply(vantVeer,1,cor,y = phenodata$metastases))
coreSet <- rownames(vantVeer)[tail(order(stat), 25)]

resOverrep <- geneSetAnalysis(
	dat = vantVeer,
	geneSets = pathways[1:4],
	analysis = analysis.customOverrepresentation(),
	coreSet = coreSet,
	adjustmentMethod = "fdr",
	signLevel=0.1)

summary(resOverrep)


###################################################
### code chunk number 14: giant_package_vignette.Snw:194-197
###################################################
pdf("overrepresentation.pdf")
plotOverrepresentation(resOverrep, subset = 1:4, aggregate = TRUE)
dev.off()


###################################################
### code chunk number 15: giant_package_vignette.Snw:200-201 (eval = FALSE)
###################################################
## plotOverrepresentation(resOverrep, aggregate = TRUE)


###################################################
### code chunk number 16: giant_package_vignette.Snw:215-229 (eval = FALSE)
###################################################
## resUncertainty <- evaluateGeneSetUncertainty(
## 	#parameters in ...
## 	labs = phenodata$metastases,
## 	numSamples = 1000,
## 	#parameters for evaluateGeneSetUncertainty
## 	dat = vantVeer,
## 	geneSet = pathways[[3]],
## 	analysis = analysis.averageCorrelation(),
## 	numSamplesUncertainty = 100,
## 	k = seq(0.1,0.9,by = 0.1))
## 
## plot(resUncertainty,
## 	main = names(pathways[3]),
## 	addMinimalStability = TRUE)


###################################################
### code chunk number 17: giant_package_vignette.Snw:231-246
###################################################
resUncertainty <- evaluateGeneSetUncertainty(
	#parameters in ...
	labs = phenodata$metastases,
	numSamples = 1000,
	#parameters for evaluateGeneSetUncertainty
	dat = vantVeer,
	geneSet = pathways[[3]],
	analysis = analysis.averageCorrelation(),
	numSamplesUncertainty = 100,
	k = seq(0.1,0.9,by = 0.1))
pdf("uncertainty.pdf")
plot(resUncertainty,
	main = names(pathways[3]),
	addMinimalStability = TRUE)
dev.off()


###################################################
### code chunk number 18: giant_package_vignette.Snw:265-270
###################################################
myGLS <- function(dat, labs, method = "pearson"){
	return(apply(dat, 1, function(x){
			cor(x = x, y = labs, method = method)
		}))
}


###################################################
### code chunk number 19: giant_package_vignette.Snw:277-280
###################################################
myGSS <- function(x, geneSetIndices){
    return(mean(x[geneSetIndices]))
}


###################################################
### code chunk number 20: giant_package_vignette.Snw:288-302
###################################################
myAnalysis <- function(){
	return(gsAnalysis(name = "myAnalysis",
		gls = "myGLS", 
		glsParameterNames = c("labs", "method"),
		transformation = "abs", 
		transformationParameterNames = NULL,
		gss = "myGSS", 
		gssParameterNames = NULL,
		globalStat = NULL,
		globalStatParameterNames = NULL, 
		significance = "significance.sampling",
		significanceParameterNames = c("numSamples"),
		testAlternative = "greater"))
}


###################################################
### code chunk number 21: giant_package_vignette.Snw:310-320
###################################################
myResult <- geneSetAnalysis(
	labs = phenodata$metastases,
	method = "pearson",
	numSamples = 100,
	dat = vantVeer,
	geneSets = pathways,
	analysis = myAnalysis(),
	adjustmentMethod = "fdr")

hist(myResult)


