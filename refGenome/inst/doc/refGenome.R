### R code from vignette source 'refGenome.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: refGenome.Rnw:36-37
###################################################
options(width=60)


###################################################
### code chunk number 2: refGenome.Rnw:80-83
###################################################
library(refGenome)
beg <- ensemblGenome()
basedir(beg) <- system.file("extdata", package="refGenome")


###################################################
### code chunk number 3: refGenome.Rnw:138-141
###################################################
ens_gtf <- "hs.ensembl.62.small.gtf"
read.gtf(beg, ens_gtf)
beg


###################################################
### code chunk number 4: refGenome.Rnw:198-204 (eval = FALSE)
###################################################
## uc <- ucscGenome()
## basedir(uc) <- "/my/ucsc/basedir"
## read.gtf(uc, "ucsc_knownGene.gtf")
## addXref(uc, "kgXref.csv")
## addEnsembl(uc, "knownToEnsembl.csv")
## addIsoforms(uc, "ucsc_knownisoforms.csv")


###################################################
### code chunk number 5: refGenome.Rnw:216-220
###################################################
ucfile <- system.file("extdata", "hs.ucsc.small.RData", package="refGenome")
uc <- loadGenome(ucfile)
ensfile <- system.file("extdata", "hs.ensembl.62.small.RData", package="refGenome")
ens <- loadGenome(ensfile)


###################################################
### code chunk number 6: refGenome.Rnw:235-236
###################################################
tableSeqids(ens)


###################################################
### code chunk number 7: refGenome.Rnw:243-245
###################################################
en1 <- extractSeqids(ens, "^1$")
en1


###################################################
### code chunk number 8: refGenome.Rnw:258-260
###################################################
ensPrimAssembly()
ucPrimAssembly()


###################################################
### code chunk number 9: refGenome.Rnw:264-268
###################################################
enpa<-extractSeqids(ens,ensPrimAssembly())
tableSeqids(enpa)
ucpa<-extractSeqids(uc,ucPrimAssembly())
tableSeqids(ucpa)


###################################################
### code chunk number 10: refGenome.Rnw:275-278
###################################################
tableFeatures(enpa)
enpf<-extractFeature(enpa,"exon")
enpf


###################################################
### code chunk number 11: refGenome.Rnw:291-293
###################################################
dxe <- extractByGeneName(enpa, "DDX11L1")
dxu <- extractByGeneName(ucpa, "DDX11L1")


###################################################
### code chunk number 12: refGenome.Rnw:302-304
###################################################
dxe <- extractByGeneId(enpa, "ENSG00000223972")
dxu <- extractByGeneId(ucpa, "ENSG00000223972")


###################################################
### code chunk number 13: refGenome.Rnw:309-311
###################################################
tableTranscript.id(enpa)
tableTranscript.id(ucpa)


###################################################
### code chunk number 14: refGenome.Rnw:315-317
###################################################
extractTranscript(ens, "ENST00000456328")
extractTranscript(uc, "uc010nxr.1")


###################################################
### code chunk number 15: refGenome.Rnw:334-338
###################################################
gpe <- getGenePositions(ens)
gpe
gpu <- getGenePositions(uc)
gpu


###################################################
### code chunk number 16: refGenome.Rnw:361-363
###################################################
enex <- refExons(ens)
ucex <- refExons(uc)


###################################################
### code chunk number 17: refGenome.Rnw:366-367
###################################################
enex


###################################################
### code chunk number 18: refGenome.Rnw:388-392
###################################################
jens <- getSpliceTable(ens)
jens
juc <- getSpliceTable(uc)
juc


###################################################
### code chunk number 19: refGenome.Rnw:411-417
###################################################
ujens <- unifyJuncs(jens)
ujuc <- unifyJuncs(juc)
jeg <- getGenePositions(jens)
jug <- getGenePositions(juc)
ujens
jug


###################################################
### code chunk number 20: refGenome.Rnw:447-456
###################################################
qry<-data.frame(
                  id=1:6,
                  start=c(10,18,61,78,82,110),
                  end=c(15,22,63,87,90,120))
ref<-data.frame(
                  id=1:5,
                  start=c(20,40,60,80,100),
                  end=c(25,45,65,85,105))
overlap(qry,ref)


###################################################
### code chunk number 21: refGenome.Rnw:531-541
###################################################
# + + + + + + + + + + + + + + + + + + #
# A) Example query  data
# + + + + + + + + + + + + + + + + + + #
##                          1       2       3       4       5       6       7 ##
qry <- data.frame(id = 1:7, seqid = "1",
            lstart = c(10100L, 11800L, 12220L, 12220L, 12220L, 32000L, 40000L),
            lend =   c(10100L, 12000L, 12225L, 12227L, 12227L, 32100L, 40100L),
            rstart = c(10200L, 12200L, 12057L, 12613L, 12650L, 32200L, 40200L),
            rend =   c(10300L, 12250L, 12179L, 12620L, 12700L, 32300L, 40300L))
##                          1       2       3       4       5       6       7 ##


###################################################
### code chunk number 22: refGenome.Rnw:550-556
###################################################
ensfile <- system.file("extdata", "hs.ensembl.62.small.RData",
                                            package="refGenome")
# Load Ensembl genome
ens <- loadGenome(ensfile)
# Calculate junction positions:
junc <- getSpliceTable(ens)


###################################################
### code chunk number 23: refGenome.Rnw:563-564
###################################################
res <- overlapJuncs(qry, junc)


###################################################
### code chunk number 24: refGenome.Rnw:596-629 (eval = FALSE)
###################################################
## library(refGenome)
## endir <- "/.../refGenomes/hsEns76"
## 
## # + + + + + + + + + + + + + + + + + + + + + + + + #
## # Read GTF
## # + + + + + + + + + + + + + + + + + + + + + + + + #
## en76 <- ensemblGenome()
## basedir(en76) <- endir
## read.gtf(en76, "Homo_sapiens.GRCh38.76.gtf")
## saveGenome(en76, "en76.RData")
## 
## # + + + + + + + + + + + + + + + + + + + + + + + + #
## # Extract primary assembly
## # + + + + + + + + + + + + + + + + + + + + + + + + #
## enpa76 <- extractSeqids(en76, ensPrimAssembly())
## saveGenome(enpa76, "enpa76.RData")
## 
## # + + + + + + + + + + + + + + + + + + + + + + + + #
## # Extract Exons
## # + + + + + + + + + + + + + + + + + + + + + + + + #
## enex76 <- refExons(enpa76)
## saveGenome(enex76, "enex76.RData")
## 
## # + + + + + + + + + + + + + + + + + + + + + + + + #
## # Extract Junctions
## # + + + + + + + + + + + + + + + + + + + + + + + + #
## enjc76 <- getSpliceTable(enpa76)
## saveGenome(enjc76, "enjc76.RData")
## 
## # + + + + + + + + + + + + + + + + + + + + + + + + #
## # Extract data.frame
## # + + + + + + + + + + + + + + + + + + + + + + + + #
## # enju76 <- unifyJuncs(enjc76)


