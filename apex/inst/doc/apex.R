## ----setup, echo=FALSE---------------------------------------------------
# set global chunk options: images will be 7x5 inches
knitr::opts_chunk$set(fig.width=7, fig.height=10, fig.path="figs/")
options(digits = 4)

## ----install, eval=FALSE-------------------------------------------------
#  library(devtools)
#  install_github("thibautjombart/apex")

## ----install2, eval=FALSE------------------------------------------------
#  install.packages("apex")

## ----load----------------------------------------------------------------
library("apex")

## ----readfiles-----------------------------------------------------------
## get address of the file within apex
files <- dir(system.file(package="apex"),patter="patr", full=TRUE)
files # this will change on your computer

## read these files
x <- read.multiFASTA(files)
x
names(x@dna) # names of the genes
par(mar=c(6,11,2,1))
plot(x)

## ----readfiles phyDat----------------------------------------------------
z <- read.multiphyDat(files, format="fasta")
z

## ----multidnaDef---------------------------------------------------------
getClassDef("multidna")

## ----multidnaclass-------------------------------------------------------
## empty object
new("multidna")

## using a list of genes as input
data(woodmouse)
woodmouse
genes <- list(gene1=woodmouse[,1:500], gene2=woodmouse[,501:965])
x <- new("multidna", genes)
x

## access the various slots
getNumInd(x) # The number of individuals
getNumLoci(x) # The number of loci
getLocusNames(x) # The names of the loci
getSequenceNames(x) # A list of the names of the sequences at each locus

getSequences(x) # A list of all loci
getSequences(x, loci = 2, simplify = FALSE) # Just the second locus (a single element list)
getSequences(x, loci = "gene1") # Just the first locus as a DNAbin object

## compare the input dataset and the new multidna
par(mfrow=c(3,1), mar=c(6,6,2,1))
image(woodmouse)
image(as.matrix(getSequences(x, 1)))
image(as.matrix(getSequences(x, 2)))

## same but with missing sequences and wrong order
genes <- list(gene1=woodmouse[,1:500], gene2=woodmouse[c(5:1,14:15),501:965])
x <- new("multidna", genes)
x
par(mar=c(6,6,2,1))
plot(x)

## ----multiphyDatDef------------------------------------------------------
getClassDef("multiphyDat")

## ----multiphyDatclass----------------------------------------------------
data(Laurasiatherian)
Laurasiatherian

## empty object
new("multiphyDat")

## simple conversion after artificially splitting data into 2 genes
genes <- list(gene1=subset(Laurasiatherian,,1:1600, FALSE),
      	 gene2=subset(Laurasiatherian,,1601:3179, FALSE))
x <- new("multiphyDat", genes, type="DNA")
x

## ----handling------------------------------------------------------------
files <- dir(system.file(package="apex"),patter="patr", full=TRUE)
files

## read files
x <- read.multiFASTA(files)
x
par(mar=c(6,11,2,1))
plot(x)

## subset
plot(x[1:3,2:4])

## ----concat, fig.width=12, fig.height=7----------------------------------
## concatenate
y <- concatenate(x)
y
par(mar=c(5,8,2,1))
image(y)

## concatenate multiphyDat object
z <- multidna2multiphyDat(x)
u <- concatenate(z)
u
tree <- pratchet(u, trace=0)
plot(tree, "u")

## ----gettree-------------------------------------------------------------
## make trees, default parameters
trees <- getTree(x)
trees

## ----hidePlotMultiPhylo, echo=TRUE,eval=FALSE----------------------------
#  plot(trees, 4, type="unrooted")

## ----plotMultiPhylo, echo=FALSE,eval=TRUE--------------------------------
par(mfrow=c(2,2)); for(i in 1:length(trees))plot(trees[[i]], type="unr")

## ----plotPhyloSingle, echo=FALSE,eval=TRUE-------------------------------
## make one single tree based on concatenated genes
tree <- getTree(x, pool=TRUE)
tree
par(mfrow=c(1,1))
plot(tree, type="unrooted")

## ----pmlPart, eval=FALSE-------------------------------------------------
#  ## input object
#  z
#  ## build trees
#  pp <- pmlPart(bf ~ edge + nni, z, control = pml.control(trace = 0))
#  pp
#  ## convert trees for plotting
#  trees <- pmlPart2multiPhylo(pp)

## ----hidePlotMultiPhylo2, echo=TRUE,eval=FALSE---------------------------
#  plot(trees, 4)

## ----plotPmlPart, echo=FALSE,eval=FALSE----------------------------------
#  par(mfrow=c(2,2)); for(i in 1:length(trees))plot(trees[[i]])

## ----export--------------------------------------------------------------
## find source files in apex
library(adegenet)
files <- dir(system.file(package="apex"),patter="patr", full=TRUE)

## import data
x <- read.multiFASTA(files)
x

## extract SNPs and export to genind
obj1 <- multidna2genind(x)
obj1


## ----mlst----------------------------------------------------------------
obj3 <- multidna2genind(x, mlst=TRUE)
obj3

## alleles of the first locus (=sequences)
alleles(obj3)[[1]]

