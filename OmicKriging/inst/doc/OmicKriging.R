### R code from vignette source 'OmicKriging.Rnw'

###################################################
### code chunk number 1: OmicKriging.Rnw:43-52
###################################################
library(OmicKriging)

binaryFile <- system.file(package = "OmicKriging",
                "doc/vignette_data/ig_genotypes.grm.bin")
binaryFileBase <- substr(binaryFile,1, nchar(binaryFile) - 4)
expressionFile <- system.file(package = "OmicKriging",
                    "doc/vignette_data/ig_gene_subset.txt.gz") 
phenotypeFile <- system.file(package = "OmicKriging",
                  "doc/vignette_data/ig_pheno.txt") 


###################################################
### code chunk number 2: OmicKriging.Rnw:58-59
###################################################
pheno <- read.table(phenotypeFile, header = T)


###################################################
### code chunk number 3: OmicKriging.Rnw:65-66
###################################################
grmMat <- read_GRMBin(binaryFileBase)


###################################################
### code chunk number 4: OmicKriging.Rnw:75-76
###################################################
gxmMat <- make_GXM(expFile = expressionFile)


###################################################
### code chunk number 5: OmicKriging.Rnw:86-91
###################################################
pcMatXM <- make_PCs_irlba(gxmMat, n.top = 10)

pcMatGM <- make_PCs_irlba(grmMat, n.top = 10)

pcMat <- cbind(pcMatGM, pcMatXM[match(rownames(pcMatGM), rownames(pcMatXM)),])


###################################################
### code chunk number 6: OmicKriging.Rnw:105-111
###################################################
result <- krigr_cross_validation(pheno.df = pheno,
	cor.list = list(grmMat, gxmMat),
	h2.vec = c(0.5, 0.5),
	covar.mat = pcMat,
	ncore = 2,
	nfold = "LOOCV")


###################################################
### code chunk number 7: closeConnetions
###################################################
allCon <- showConnections()
socketCon <- as.integer(rownames(allCon)[allCon[, "class"] == "sockconn"])
sapply(socketCon, function(ii) close.connection(getConnection(ii)) )


