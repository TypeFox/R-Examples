### R code from vignette source 'polysattutorial.Rnw'

# Obtaining and Installing polysat
install.packages("combinat")
install.packages("polysat")
library(polysat)

# Getting started: at tutorial
## Creating a dataset
getwd()

simgen <- read.GeneMapper("GeneMapperExample.txt")

summary(simgen)
Samples(simgen)
Loci(simgen)
viewGenotypes(simgen, samples=paste("A", 1:20, sep=""), loci="loc1")
find.missing.gen(simgen)

Description(simgen) <- "Dataset for the tutorial"
PopNames(simgen) <- c("PopA", "PopB", "PopC")
PopInfo(simgen) <- rep(1:3, each = 100)
Usatnts(simgen) <- c(2, 3, 2)

rep(1:3, each = 100)
PopInfo(simgen)

Samples(simgen, populations = "PopA")

Usatnts(simgen)

simgen <- editGenotypes(simgen, maxalleles = 4)

simgen <- estimatePloidy(simgen)

summary(simgen)

save(simgen, file="simgen.RData")

## Data analysis and export
### Genetic distances between individuals
testmat <- meandistance.matrix(simgen)
pca <- cmdscale(testmat)
mycol <- c("red", "green", "blue")
plot(pca[,1], pca[,2], col=mycol[PopInfo(simgen)],
     main = "PCA with Bruvo distance")

testmat2 <- meandistance.matrix(simgen, distmetric=Lynch.distance)

pca2 <- cmdscale(testmat2)
plot(pca2[,1], pca2[,2], col=rep(c("red", "green", "blue"), each=100),
     main = "PCA with Lynch distance")

### Working with subsets of data
simgen2 <- deleteSamples(simgen, c("B59", "C30"))
simgen2 <- deleteLoci(simgen2, "loc2")
summary(simgen2)

samToUse <- Samples(simgen2, populations=c("PopA", "PopB"), ploidies=4)
exclude <- c("A50", "A78", "B25", "B60", "B81")
samToUse <- samToUse[!samToUse %in% exclude]
samToUse

summary(simgen2[samToUse, "loc1"])

testmat3 <- meandistance.matrix(simgen2, samples = samToUse,
                                distmetric = Lynch.distance,
                                progress= FALSE)
pca3 <- cmdscale(testmat3)
plot(pca3[,1], pca3[,2], col=c("red", "blue")[PopInfo(simgen2)[samToUse]])


### Population statistics
simfreq <- deSilvaFreq(simgen, self = 0.1, initNull = 0.01,
                       samples = Samples(simgen, ploidies = 4))

simfreq
simFst <- calcFst(simfreq)
simFst
simFst12 <- calcFst(simfreq, loci=c("loc1", "loc2"))
simFst12

### Genotype data export
write.Structure(simgen, ploidy = 4, file="simgenStruct.txt")

# How data are stored in polysat
## The "genambig" class
showClass("genambig")

mysamples <- c("indA", "indB", "indC", "indD", "indE", "indF")
myloci <- c("loc1", "loc2", "loc3")
mydataset <- new("genambig", samples=mysamples, loci=myloci)

mydataset

?Samples

Loci(mydataset)
Loci(mydataset) <- c("L1", "L2", "L3")
Loci(mydataset)
Samples(mydataset)
Samples(mydataset)[3] <- "indC1"
Samples(mydataset)
PopNames(mydataset) <- c("Yosemite", "Sequoia")
PopInfo(mydataset) <- c(1,1,1,2,2,2)
PopInfo(mydataset)
PopNum(mydataset, "Yosemite")
PopNum(mydataset, "Sequoia") <- 3
PopNames(mydataset)
PopInfo(mydataset)
Ploidies(mydataset) <- c(4,4,4,4,4,6)
Ploidies(mydataset)

Ploidies(mydataset)["indC1",] <- 6
Ploidies(mydataset)
Usatnts(mydataset) <- c(2,2,2)
Usatnts(mydataset)
Description(mydataset) <- "Tutorial, part 2."
Description(mydataset)
Genotypes(mydataset, loci="L1") <- list(c(122, 124, 128), c(124,126),
                     c(120,126,128,130), c(122,124,130), c(128,130,132),
                     c(126,130))
Genotype(mydataset, "indB", "L3") <- c(150, 154, 160)
Genotypes(mydataset)
Genotype(mydataset, "indD", "L1")
Missing(mydataset)
Missing(mydataset) <- -1
Genotypes(mydataset)

mydataset@Genotypes
mydataset@Genotypes[["indB","L1"]]

isMissing(mydataset, "indA", "L2")
isMissing(mydataset, "indA", "L1")
isMissing(mydataset)

moredata <- new("genambig", samples=c("indG", "indH"), loci=Loci(mydataset))
Usatnts(moredata) <- Usatnts(mydataset)
Description(moredata) <- Description(mydataset)
PopNames(moredata) <- "Kings Canyon"
PopInfo(moredata) <- c(1,1)
Ploidies(moredata) <- c(4,4)
Missing(moredata) <- Missing(mydataset)
Genotypes(moredata, loci="L1") <- list(c(126,130,136,138), c(124,126,128))
mydataset2 <- merge(mydataset, moredata)
mydataset2

## How ploidy is stored: "ploidysuper" and subclasses
ploidyexample <- new("genambig")
Samples(ploidyexample)
Loci(ploidyexample)
Ploidies(ploidyexample)
ploidyexample <- reformatPloidies(ploidyexample, output="locus")
Ploidies(ploidyexample)
ploidyexample <- reformatPloidies(ploidyexample, output="sample")
Ploidies(ploidyexample)
ploidyexample <- reformatPloidies(ploidyexample, output="one")
Ploidies(ploidyexample)
Ploidies(ploidyexample) <- 4
ploidyexample <- reformatPloidies(ploidyexample, output="matrix")
Ploidies(ploidyexample)

Ploidies(ploidyexample)["ind1", "loc1"]

Ploidies(ploidyexample, "ind1", "loc1")
ploidyexample <- reformatPloidies(ploidyexample, output="one")
Ploidies(ploidyexample, "ind1", "loc1")

## The "gendata" and "genbinary" classes
simgenB <- genambig.to.genbinary(simgen)
Genotypes(simgenB, samples=paste("A", 1:20, sep=""), loci="loc1")
PopInfo(simgenB)[Samples(simgenB, ploidies=2)]

write.table(Genotypes(simgenB), file="simBinaryData.txt")

# Functions for autopolyploid data
## Data import
GDdata <- read.GenoDive("genodiveExample.txt")
Structdata <- read.Structure("structureExample.txt", ploidy = 8)
Spagdata <- read.SPAGeDi("spagediExample.txt")
PDdata <- read.POPDIST(c("POPDISTexample1.txt", "POPDISTexample2.txt"))

GMdata <- read.GeneMapper(c("GeneMapperCBA15.txt",
                            "GeneMapperCBA23.txt",
                            "GeneMapperCBA28.txt"))

# view the format
read.table("STRandExample.txt", sep="\t", header=TRUE)
# import the data
STRdata <- read.STRand("STRandExample.txt")
Samples(STRdata)
PopNames(STRdata)

domdata <- read.table("dominantExample.txt", header=TRUE,
                      sep="\t", row.names=1)

domdata
domdata <- as.matrix(domdata)
PAdata <- new("genbinary", samples=c("ind1", "ind2", "ind3"),
              loci=c("ABC1", "ABC2"))
Genotypes(PAdata) <- domdata

PopInfo(PAdata) <- c(1,1,2)
PAdata <- genbinary.to.genambig(PAdata)

## Data export
myexcol <- array(c(rep(0:1, each=150), seq(0.1, 30, by=0.1)), dim=c(300,2),
                 dimnames = list(Samples(simgen), c("PopFlag", "Something")))
myexcol[1:10,]
write.Structure(simgen, ploidy=4, file="simgenStruct2.txt",
                writepopinfo = FALSE, extracols = myexcol,
                missingout = -1)

write.GenoDive(simgen, file="simgenGD.txt")

write.SPAGeDi(simgen, file="simgenSpag.txt")

write.POPDIST(simgen, samples = Samples(simgen, ploidies=4),
              file = "simgenPOPDIST.txt")

write.GeneMapper(simgen, file="simgenGM.txt")

simgenPA <- genambig.to.genbinary(simgen)
write.table(Genotypes(simgenPA), file="simgenPA.txt", quote=FALSE,
            sep = ",")


## Individual-level statistics
### Estimating and exporting ploidies
write.table(data.frame(Ploidies(simgen), row.names=Samples(simgen)),
            file="simgenPloidies.txt")

### Inter-individual distances
testmat4 <- meandistance.matrix2(simgen, samples=samToUse, freq=simfreq,
                                 self=0.2)
pca4 <- cmdscale(testmat4)
plot(pca4[,1], pca4[,2], col=c("red", "blue")[PopInfo(simgen)[samToUse]],
     main="Bruvo distance with meandistance.matrix2")

hist(as.vector(testmat))

hist(as.vector(testmat2))

write.table(testmat2, file="simgenDistMat.txt")

subsamples <- Samples(simgen, populations=1)
subsamples <- subsamples[!isMissing(simgen, subsamples, "loc1") &
                         !isMissing(simgen, subsamples, "loc2") &
                         !isMissing(simgen, subsamples, "loc3")]
Larray <- meandistance.matrix(simgen, samples=subsamples,
                              progress=FALSE,
                 distmetric=Lynch.distance, all.distances=TRUE)[[1]]
mdist1.2 <- meandist.from.array(Larray, loci=c("loc1","loc2"))
mdist2.3 <- meandist.from.array(Larray, loci=c("loc2","loc3"))
mdist1.3 <- meandist.from.array(Larray, loci=c("loc1","loc3"))

### Determining groups of asexually-related samples
clones <- assignClones(testmat, samples=paste("A", 1:100, sep=""),
                       threshold=0.2)
clones

## Population statistics
### Allele diversity and frequencies
simal <- alleleDiversity(simgen)
simal$counts
simal$alleles[["PopA","loc1"]]

simFst
simfreqSimple <- simpleFreq(simgen, samples = Samples(simgen, ploidies=4))
simFstSimple <- calcFst(simfreqSimple)
simFstSimple

write.freq.SPAGeDi(simfreq, usatnts=Usatnts(simgen), file="SPAGfreq.txt")

gpsimfreq <- freq.to.genpop(simfreq)

### Genotype frequencies
testmat5 <- meandistance.matrix(simgen, all.distances=TRUE)
simdiv <- genotypeDiversity(simgen, d=testmat5, threshold=0.2, index=Shannon)
simdiv
simdiv2 <- genotypeDiversity(simgen, d=testmat5, threshold=0.2, index=Simpson)
simdiv2

simdiv2var <- genotypeDiversity(simgen, d=testmat5, threshold=0.2,
                                index=Simpson.var)
simdiv2 - 2*simdiv2var^0.5
simdiv2 + 2*simdiv2var^0.5

# Functions for allopolyploid data
## Data import and export
ATdata <- read.ATetra("ATetraExample.txt")
Tetdata <- read.Tetrasat("tetrasatExample.txt")

write.ATetra(simgen, samples=Samples(simgen, ploidies=4), file="simgenAT.txt")
write.Tetrasat(simgen, samples=Samples(simgen, ploidies=4),
               file="simgenTet.txt")

# Treating microsatellite alleles as dominant markers
Present(simgenB) <- "P"
Absent(simgenB) <- 2
Missing(simgenB) <- 0
Genotypes(simgenB)[1:10, 1:6]

genmat <- Genotypes(simgenB)
dimnames(genmat)[[2]] <- paste("M", 1:dim(genmat)[2], sep="")
genmat[1:10, 1:10]


