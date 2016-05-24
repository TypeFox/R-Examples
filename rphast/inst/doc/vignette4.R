### R code from vignette source 'vignette4.Rnw'

###################################################
### code chunk number 1: readHar1
###################################################
getOption("SweaveHooks")[["fig"]]()
require("rphast")
exampleArchive <- system.file("extdata", "examples.zip", package="rphast")
unzip(exampleArchive, c("HAR_001_neutral_SSREV.mod", "HAR_001.fa"))
neutralMod <- read.tm("HAR_001_neutral_SSREV.mod")
align <- read.msa("HAR_001.fa")


###################################################
### code chunk number 2: HAReasy
###################################################
getOption("SweaveHooks")[["fig"]]()
nuc.results <- bgc.nucleotide.tests(align, neutralMod, "hg18")
nuc.results


###################################################
### code chunk number 3: countMuts
###################################################
getOption("SweaveHooks")[["fig"]]()
classify.muts.bgc(align, neutralMod, branch="hg18")


###################################################
### code chunk number 4: ncol
###################################################
getOption("SweaveHooks")[["fig"]]()
ncol.msa(align)


###################################################
### code chunk number 5: HAReasyBounds
###################################################
getOption("SweaveHooks")[["fig"]]()
bgc.nucleotide.tests(align, neutralMod, "hg18", bgc.limits=c(0, 2000), sel.limits=c(-2000,2000))


###################################################
### code chunk number 6: HAR
###################################################
getOption("SweaveHooks")[["fig"]]()
# note this is not portable
nullLike <- nuc.results["null", "likelihood"]
selLike <- nuc.results["sel", "likelihood"]
bgcLike <- nuc.results["bgc", "likelihood"]
bgcSelLike <- nuc.results["sel+bgc", "likelihood"]

lrs <- c(selLike - nullLike,
         bgcLike - nullLike,
         selLike - bgcLike,
         bgcLike - selLike,
         bgcSelLike - bgcLike,
         bgcSelLike - selLike)
lrs


###################################################
### code chunk number 7: HMMReadData
###################################################
getOption("SweaveHooks")[["fig"]]()
unzip(exampleArchive, c("chr1_100k_4way.maf", "placentalMammals.mod"))
align <- read.msa("chr1_100k_4way.maf")
mod <- read.tm("placentalMammals.mod")
mod$tree <- prune.tree(mod$tree, names(align), all.but=TRUE)


###################################################
### code chunk number 8: runBgcHMM
###################################################
getOption("SweaveHooks")[["fig"]]()
hmm.results <- phastBias(align, mod, foreground="hg18")
hmm.results


###################################################
### code chunk number 9: showTract
###################################################
getOption("SweaveHooks")[["fig"]]()
hmm.results$tracts
coverage.feat(hmm.results$tracts)


###################################################
### code chunk number 10: extractTract
###################################################
getOption("SweaveHooks")[["fig"]]()
tractAlign <- split.by.feature.msa(align, hmm.results$tracts)[[1]]
classify.muts.bgc(tractAlign, mod, "hg18")
classify.muts.bgc(align, mod, "hg18")


###################################################
### code chunk number 11: postPlot
###################################################
getOption("SweaveHooks")[["fig"]]()

tracks <- list(as.track.feat(hmm.results$tracts, name="gBGC tracts"),
               as.track.wig(score=(hmm.results$post.prob[,"gBGC.neutral"] + 
                                   hmm.results$post.prob[,"gBGC.conserved"]),
                            name="gBGC posterior", 
                            coord=hmm.results$post.prob$coord),
               as.track.feat(hmm.results$informative, name="informative"))
plot.track(tracks)



###################################################
### code chunk number 12: plotMsa
###################################################
getOption("SweaveHooks")[["fig"]]()
plot.msa(tractAlign, xlim=c(1083115, 1083175), pretty=TRUE)


###################################################
### code chunk number 13: reloadHAR1
###################################################
getOption("SweaveHooks")[["fig"]]()
align <- read.msa("HAR_001.fa")
neutralMod <- read.tm("HAR_001_neutral_SSREV.mod")
neutralMod


###################################################
### code chunk number 14: addSelHar1
###################################################
getOption("SweaveHooks")[["fig"]]()
neutralMod$selection <- 0


###################################################
### code chunk number 15: estimateNull
###################################################
getOption("SweaveHooks")[["fig"]]()
nullMod <- phyloFit(align, init.mod=neutralMod, no.opt=c("backgd", "branches", "ratematrix"))
nullMod


###################################################
### code chunk number 16: rescaleNull
###################################################
getOption("SweaveHooks")[["fig"]]()
phyloFit(align, init.mod=neutralMod, scale.only=TRUE, no.opt=c("backgd", "ratematrix", "sel"))


###################################################
### code chunk number 17: checkRescale
###################################################
getOption("SweaveHooks")[["fig"]]()
rescale.tree(neutralMod$tree, nullMod$selection/(1-exp(-nullMod$selection)))


###################################################
### code chunk number 18: nullLike
###################################################
getOption("SweaveHooks")[["fig"]]()
nullMod$likelihood
likelihood.msa(align, nullMod)


###################################################
### code chunk number 19: nullSelection
###################################################
getOption("SweaveHooks")[["fig"]]()
nullMod$selection


###################################################
### code chunk number 20: addAltModSel
###################################################
getOption("SweaveHooks")[["fig"]]()
initSelMod <- add.ls.mod(neutralMod, "hg18", separate.params="sel")


###################################################
### code chunk number 21: addAltModSel2
###################################################
getOption("SweaveHooks")[["fig"]]()
initSelMod2 <- add.ls.mod(neutralMod, "hg18", separate.params="sel", selection=2)


###################################################
### code chunk number 22: estimateSel
###################################################
getOption("SweaveHooks")[["fig"]]()
selMod <- phyloFit(align, init.mod=initSelMod, no.opt=c("backgd", "branches", "ratematrix"))

# print the likelihood:
selMod$likelihood
# print the global selection parameter:
selMod$selection
# print the human selection parameter:
selMod$ls.mod$selection


###################################################
### code chunk number 23: checkSel
###################################################
getOption("SweaveHooks")[["fig"]]()
apply.bgc.sel(neutralMod$rate.matrix, sel=(selMod$selection+selMod$ls.mod$selection)) - 
  selMod$ls.model$rate.matrix


###################################################
### code chunk number 24: addAltModBgc
###################################################
getOption("SweaveHooks")[["fig"]]()
initBgcMod <- add.ls.mod(neutralMod, "hg18", separate.params="bgc[0,2000]")


###################################################
### code chunk number 25: estimateBgc
###################################################
getOption("SweaveHooks")[["fig"]]()
bgcMod <- phyloFit(align, init.mod=initBgcMod, no.opt=c("backgd", "branches", "ratematrix"))

#print the likelihood, global selection, and bgc parameters:
bgcMod$likelihood
bgcMod$selection
bgcMod$ls.mod$bgc


###################################################
### code chunk number 26: addAltModBgcSel
###################################################
getOption("SweaveHooks")[["fig"]]()
initBgcSelMod <- add.ls.mod(neutralMod, "hg18", separate.params=c("bgc[0,2000]", "sel[-1000,1000]"))
bgcSelMod <- phyloFit(align, init.mod=initBgcSelMod, no.opt=c("backgd", "branches", "ratematrix"))

# print likelihood, global selection, lineage-specific selection, and bgc parameters:
bgcSelMod$likelihood
bgcSelMod$selection
bgcSelMod$alt.mod$selection
bgcSelMod$alt.mod$bgc



