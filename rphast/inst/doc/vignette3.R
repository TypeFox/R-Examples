### R code from vignette source 'vignette3.Rnw'

###################################################
### code chunk number 1: vignette3.Rnw:21-23
###################################################
options(SweaveHooks=list(fig=function()
          par(mar=c(3.1,2.1,4.1,2.1))))


###################################################
### code chunk number 2: phyloHmm1
###################################################
require("rphast")
seqnames <- c("hg18", "panTro2", "ponAbe2", "rheMac2", "equCab2", 
              "canFam2", "dasNov2", "mm9", "rn4", "cavPor3") 
exampleArchive <- system.file("extdata", "examples.zip", package="rphast")
unzip(exampleArchive, "placentalMammals.mod")
neutralMod <- read.tm("placentalMammals.mod")
neutralMod$tree <- prune.tree(neutralMod$tree, seqs=seqnames, all.but=TRUE)


###################################################
### code chunk number 3: phyloHmm2
###################################################
unzip(system.file("extdata", "NRSF.zip", package="rphast"))
mafFiles <- list.files("NRSF", pattern="*.maf", full.names=TRUE)
nrsfNames <- sub(".maf$", "", basename(mafFiles))  # remove dir and .maf 
nrsfSites <- read.feat("NRSF/NRSF.gff")
msaList <- list()
for (i in 1:length(mafFiles)) {
  smallMsa <- read.msa(mafFiles[i], seqnames=seqnames)
  smallMsa <- strip.gaps.msa(smallMsa)
  smallMsa$offset <- 0 
  feat <- nrsfSites[which(nrsfSites$feature == nrsfNames[i]),]
  if (feat$strand == "-") 
    smallMsa <- reverse.complement.msa(smallMsa)
  msaList[[nrsfNames[i]]] <- smallMsa
}
aggMsa <- concat.msa(msaList)
motifLen <- unique(sapply(msaList, ncol))
if (length(motifLen) != 1L) warning("all motifs should have same length!\n")


###################################################
### code chunk number 4: phyloHmm3
###################################################
feats <- feat(seqname="hg18", src=as.character(sapply(nrsfNames, rep, motifLen)),
              feature=rep(sprintf("site.%i", 1:motifLen), length(mafFiles)),
              start=1:ncol(aggMsa), end=1:ncol(aggMsa))
mods <- phyloFit(aggMsa, init.mod=neutralMod, no.opt="ratematrix", 
                 features=feats,
                 scale.only=TRUE, ninf.sites=10)


###################################################
### code chunk number 5: lrByMotifPosition
###################################################
getOption("SweaveHooks")[["fig"]]()
nullLike <- numeric()
lr  <- numeric()

for (i in 1:motifLen) {
  nullLike[i] <- likelihood.msa(concat.msa(lapply(msaList, `[`, cols=i)), neutralMod)
  lr[i] <- mods[[i]]$likelihood - nullLike[i]
}
barplot(lr, names.arg=1:motifLen, ylab="likelihood ratio")


###################################################
### code chunk number 6: seqLogoPlot
###################################################
getOption("SweaveHooks")[["fig"]]()
haveSeqLogo <- require("seqLogo")
if (haveSeqLogo) {
    m <- read.table("NRSF/NRSF.mtf")
    pwm <- makePWM(t(m))
    seqLogo(pwm, xfontsize=10)
} else {
    plot(c(0),c(0), type="n", xlab="", ylab="", xaxt="n", yaxt="n")
    text(0, 0, "Need to install seqLogo to see this plot")
}


###################################################
### code chunk number 7: seqLogoPlotBackgd
###################################################
getOption("SweaveHooks")[["fig"]]()
if (haveSeqLogo) {
    m <- matrix(0, nrow=21, ncol=4)
    for (i in 1:21)
        m[i,] <- mods[[i]]$backgd
    pwm <- makePWM(t(m))
    seqLogo(pwm, xfontsize=10)
} else {
    plot(c(0),c(0), type="n", xlab="", ylab="", xaxt="n", yaxt="n")
    text(0, 0, "Need to install seqLogo to see this plot")
}


###################################################
### code chunk number 8: phyloHmm7
###################################################
mods[["neutral"]] <- neutralMod

get.trans.mat <- function(lambda, state.names, motifLen) {
  transMat <- matrix(0, nrow=length(state.names), ncol=length(state.names),
                      dimnames=list(state.names, state.names))
  transMat["neutral", "site.1"] <- lambda
  transMat["neutral", "neutral"] <- 1 - lambda
  for (i in 1:(motifLen-1)) 
    transMat[sprintf("site.%i", i), sprintf("site.%i", i+1)] <- 1
  transMat[sprintf("site.%i", motifLen), "neutral"] <- 1.0
  transMat
}

lambda <- 0.0001
transMat <- get.trans.mat(lambda, names(mods), motifLen)
nrsfHmm <- hmm(transMat)


###################################################
### code chunk number 9: phyloHmm8
###################################################
simLength <- 100000
simData <- simulate.msa(mods, simLength, hmm=nrsfHmm, get.features=TRUE)


###################################################
### code chunk number 10: phyloHmm9
###################################################
hmmScores <- score.hmm(msa=simData$msa, mod=mods, hmm=nrsfHmm, viterbi=TRUE,
                       states=sprintf("site.%i", 1:motifLen))


###################################################
### code chunk number 11: phyloHmm10
###################################################
predicted <- hmmScores$in.states
correct <- simData$feats[substr(simData$feats$feature, 1, 4)=="site",]

numSitePredicted <- coverage.feat(predicted)
numSiteCorrect <- coverage.feat(correct)
numSiteOverlap <- coverage.feat(predicted, correct)
cat(numSitePredicted, numSiteCorrect, numSiteOverlap, "\n")


###################################################
### code chunk number 12: vignette3BrowserView
###################################################
getOption("SweaveHooks")[["fig"]]()

wholeRegion <- feat("hg18", src=".", feature="all",start=1, end=simLength)
sensitivity <- coverage.feat(correct, predicted)/coverage.feat(correct)
specificity <- coverage.feat(correct, predicted, wholeRegion, not=c(TRUE, TRUE, FALSE))/ 
  coverage.feat(correct, wholeRegion, not=c(TRUE, FALSE))
ppv <- coverage.feat(correct, predicted) /
  coverage.feat(predicted)
cat(specificity, sensitivity, ppv, "\n")

tracks <- list(as.track.feat(correct, "actual binding sites"),
               as.track.feat(predicted, "predicted binding sites", col="red"),
               as.track.wig(coord=hmmScores$post.prob.wig$coord,
                            score=hmmScores$post.prob.wig$post.prob,
                            name="hmm Posterior probilities",
                            smooth=FALSE, col="red", ylim=c(0,1)))
plot.track(tracks)


###################################################
### code chunk number 13: phyloHmm12
###################################################
unlink("NRSF", recursive=TRUE)


