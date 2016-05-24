###############################################################################
## Affymetrix example
###############################################################################

###############################################################################
## References:
## L.M. Cope1 , R.A. Irizarry, H.A. Jaffee, Z. Wu and T.P. Speed (2004):
## "A benchmark for Affymetrix GeneChip expression measures". 
## Bioinformatics 20(3): 323-331
##
## R.A. Irizarry1, Z. Wu and H.A. Jaffee (2006):
## "Comparison of Affymetrix GeneChip expression measures"
## Bioinformatics 22(7): 789-794
###############################################################################

###############################################################################
## Data
## Spike-in hgu95a data:
## http://www.biostat.jhsph.edu/~ririzarr/affycomp/spikein.tgz
##
## Spike-in hgu133a data:
## http://www.biostat.jhsph.edu/~ririzarr/affycomp/hgu133spikein.tgz
##
## Dilution data: 
## The links:
## http://qolotus02.genelogic.com/datasets.nsf/               (no response)
## and
## http://www.genelogic.com/media/studies/dilution.cfm        (not found)
## seem not to lead to the data any longer.
## An email to the support of genelogic remained unanswered ...
###############################################################################


library(affy)
library(affycomp) ## Version >= 1.19.4

###################
## replace with your path to hgu95a data!!!
PATH <- "./spikein"
###################

fn <- list.celfiles(path = PATH, full.names=TRUE)
data(spikein.phenodata)
(spikein.hgu95a <- read.affybatch(filenames = fn, 
                                  phenoData = spikein.phenodata))

###################
## replace with your path to hgu133a data!!!
PATH <- "./hgu133spikein"
###################

fn <- list.celfiles(path = PATH, full.names=TRUE)
fn <- fn[c(seq(1,40,3), seq(2, 41, 3), seq(3, 42, 3))]
fn <- fn[c(6:14, 1:5, 20:28, 15:19, 34:42, 29:33)]
data(hgu133a.spikein.phenodata)

## Attention:
## Order of filenames in fn has to be identical to 
## sampleNames(hgu133a.spikein.phenodata)!!!

(spikein.hgu133a <- read.affybatch(filenames = fn,
                                  phenoData = hgu133a.spikein.phenodata))


###########################################################
## Is normal location and scale appropriate for the log-data?
## As there is no way to test for approximate normality, we use the Kolmogorov
## distance as indicator and compute the minimum Kolmogorov distance between
## the log-data and the set of (univariate) normal distributions.
###########################################################
library(RobLoxBioC)

## takes more than 100 min on Intel P9500 (64bit Linux, 4 GByte RAM)
system.time(minKD.hgu95a <- KolmogorovMinDist(spikein.hgu95a, Norm()))
## takes more than 130 min on Intel P9500 (64bit Linux, 4 GByte RAM)
system.time(minKD.hgu133a <- KolmogorovMinDist(spikein.hgu133a, Norm()))

## load the results from R-forge ...
con <- url("http://robast.r-forge.r-project.org/data/minKD_hgu95a.RData")
load(file = con)
close(con)
con <- url("http://robast.r-forge.r-project.org/data/minKD_hgu133a.RData")
load(file = con)
close(con)

boxplot(as.data.frame(minKD.hgu95a$dist), main = "HGU95a")
boxplot(as.data.frame(minKD.hgu133a$dist), main = "HGU133a")

## Table 1 in Kohl and Deigner (2010)
table(minKD.hgu95a$n)/59
table(minKD.hgu133a$n)/42

###########################################################
## Comparison with normal (pseudo) random samples
###########################################################

## takes more than 90 min on Intel P9500 (64bit Linux, 4 GByte RAM)
ns <- c(5:20, 69)
M <- length(ns)
minKD.norm <- matrix(NA, nrow = 50000, ncol = M)
for(i in seq_len(M)){
    print(i)
    temp <- matrix(rnorm(50000*ns[i]), ncol = ns[i])
    minKD.norm[,i] <- KolmogorovMinDist(temp, Norm())$dist
}
colnames(minKD.norm) <- ns

## load the results from R-forge
con <- url("http://robast.r-forge.r-project.org/data/minKD_norm.RData")
load(file = con)
close(con)

x11(width = 14)
par(mfrow = c(1, 3))
res <- split(as.vector(minKD.hgu95a$dist), as.vector(minKD.hgu95a$n))
uni.n <- sort(unique(as.vector(minKD.hgu95a$n)))
boxplot(res, main = "HGU95a", ylim = c(0, 0.45), 
        ylab = "minimum Kolmogorov distance", xlab = "sample size")
lines(1:length(uni.n), 1/(2*uni.n), col = "orange", lwd = 2)
legend("topright", legend = "minimal possible distance", fill = "orange")

res <- split(as.vector(minKD.hgu133a$dist), as.vector(minKD.hgu133a$n))
uni.n <- sort(unique(as.vector(minKD.hgu133a$n)))
boxplot(res, main = "HGU133a", ylim = c(0, 0.45), 
        ylab = "minimum Kolmogorov distance", xlab = "sample size")
lines(1:length(uni.n), 1/(2*uni.n), col = "orange", lwd = 2)
legend("topright", legend = "minimal possible distance", fill = "orange")

boxplot(as.data.frame(minKD.norm), main = "Normal samples", ylim = c(0, 0.45), 
        ylab = "minimum Kolmogorov distance", xlab = "sample size")
lines(1:length(ns), 1/(2*ns), col = "orange", lwd = 2)
legend("topright", legend = "minimal possible distance", fill = "orange")

#######################################
## Figure 2 in Kohl and Deigner (2010)
#######################################
res1 <- split(as.vector(minKD.hgu95a$dist), as.vector(minKD.hgu95a$n))
res2 <- split(as.vector(minKD.hgu133a$dist), as.vector(minKD.hgu133a$n))
res3 <- lapply(as.data.frame(minKD.norm[,c(2:12,16,17)]), function(x) x)
uni.n <- c(as.integer(names(res1)), as.integer(names(res2)), as.integer(names(res3)))

#setEPS(height = 6, width = 9)
#postscript(file = "Figure2.eps")
par(mar = c(4, 4, 3, 1))
plot(0, 0, type = "n", ylim = c(0, 0.49), xlim = c(0.5, 16.5), 
     panel.first = abline(h = seq(0, 0.45, by = 0.05), lty = 2, col = "grey"), 
     main = "Minimum Kolmogorov distance", 
     ylab = "minimum Kolmogorov distance", 
     xlab = "sample size", axes = FALSE)
axis(1, c(1:5, 7:9, 11:16), labels = uni.n[c(8:12, 6,11,12,6,8:12)], cex.axis = 0.6)
axis(2, seq(0, 0.45, by = 0.05), labels = seq(0, 0.45, by = 0.05), las = 2,
     cex.axis = 0.8)
box()
boxplot(c(res1[8:12], res2[c(3,7,8)], res3[c(6,8:12)]), at = c(1:5, 7:9, 11:16), add = TRUE, pch = 20, 
        names = FALSE, axes = FALSE, 
        col = c(rep(NA, 3), grey(0.6), NA, grey(0.4), NA, NA, grey(0.4), rep(NA, 3),
                grey(0.6), NA))
abline(v = c(6, 10), lwd = 1.5)
text(c(3, 8, 13.5), rep(0.48, 3), labels = c("HGU95A", "HGU133A", "Normal Samples"),
     font = 2)
lines(1:5, 1/(2*uni.n[8:12]), lwd = 2)
lines(7:9, 1/(2*uni.n[c(6,8,12)]), lwd = 2)
lines(11:16, 1/(2*uni.n[c(6,8:12)]), lwd = 2)
legend("bottomleft", legend = "minimal possible distance", lty = 1, lwd = 2, cex = 0.8)
abline(h = c(0.1, 0.15), lty = 2, lwd = 1.5)
#dev.off()

## Comparison of median distances
## Table 2 in Kohl and Deigner (2010)
round(sapply(res1, quantile, prob = 0.95) - sapply(res3, quantile, prob = 0.95), 4)
round(sapply(res2, quantile, prob = 0.95) - sapply(res3, quantile, prob = 0.95)[-c(1,2,4,7)], 4)

round(sapply(res1, median) - sapply(res3, median), 4)
round(sapply(res2, median) - sapply(res3, median)[-c(1,2,4,7)], 4)


###########################################################
## assessments for MAS 5.0 and RMA including dilution data from package affycomp
###########################################################
data(mas5.assessment)
data(rma.assessment)
data(mas5.assessment.133)
data(rma.assessment.133)

## just for timing ...
## on Intel P9500 (64 bit Linux, 4 GByte RAM)

## about 510 sec
#system.time(mas5.res <- mas5(spikein.hgu95a))

## about 570 sec
#system.time(mas5.res.133 <- mas5(spikein.hgu133a))

## Implementation of function mas5 in affy package could clearly be improved ...
## by parallelizing function affy::tukey.biweight (similar to roblox) and by 
## modifying function affy::pmcorrect.mas
## function (object, contrast.tau = 0.03, scale.tau = 10, delta = 2^(-20)) 
## {
##    ...
##    for (i in 1:ncol(diff)) {
##        ...
##
##        the following line is quite slow ... much faster pmax.int!
##        !!!
##        pm.corrected <- apply(cbind(pps.pm - pps.im, delta), 1, max)
##        !!!
##        ...
##    }
##    return(diff)
## }



## about 30 sec
#system.time(rma.res <- rma(spikein.hgu95a))

## about 26 sec
#system.time(rma.res.133 <- rma(spikein.hgu133a))


###########################################################
## Example 1: Analogous to "classical" MAS 5.0
## both computations take about 55 sec on Intel P9500 (64bit Linux, 4 GByte RAM)
###########################################################
## hgu95a
system.time(eset.hgu95a <- robloxbioc(spikein.hgu95a, normalize = TRUE, add.constant = 0))
eset.hgu95a.log2 <- eset.hgu95a
exprs(eset.hgu95a.log2) <- log2(exprs(eset.hgu95a))
roblox.hgu95a <- assessSpikeIn(eset.hgu95a.log2, method.name = "roblox")
roblox.hgu95a.2 <- assessSpikeIn2(eset.hgu95a.log2, method.name = "roblox")

## hgu133a
system.time(eset.hgu133a <- robloxbioc(spikein.hgu133a, normalize = TRUE, add.constant = 0))
eset.hgu133a.log2 <- eset.hgu133a
exprs(eset.hgu133a.log2) <- log2(exprs(eset.hgu133a))
roblox.hgu133a <- assessSpikeIn(eset.hgu133a.log2, method.name = "roblox")
roblox.hgu133a.2 <- assessSpikeIn2(eset.hgu133a.log2, method.name = "roblox")


###########################################################
## Example 2: MAS 5.0 + 32
## both computations take about 55 sec on Intel P9500 (64bit Linux, 4 GByte RAM)
###########################################################
## hgu95a
system.time(eset.hgu95a32 <- robloxbioc(spikein.hgu95a, normalize = TRUE, add.constant = 32))
eset.hgu95a.log232 <- eset.hgu95a32
exprs(eset.hgu95a.log232) <- log2(exprs(eset.hgu95a32))
roblox.hgu95a32 <- assessSpikeIn(eset.hgu95a.log232, method.name = "roblox + 32")
roblox.hgu95a32.2 <- assessSpikeIn2(eset.hgu95a.log232, method.name = "roblox + 32")

## hgu133a
system.time(eset.hgu133a32 <- robloxbioc(spikein.hgu133a, normalize = TRUE, add.constant = 32))
eset.hgu133a.log232 <- eset.hgu133a32
exprs(eset.hgu133a.log232) <- log2(exprs(eset.hgu133a32))
roblox.hgu133a32 <- assessSpikeIn(eset.hgu133a.log232, method.name = "roblox + 32")
roblox.hgu133a32.2 <- assessSpikeIn2(eset.hgu133a.log232, method.name = "roblox + 32")


###########################################################
## Example 3: MAS 5.0 with PM only
## both computations take about 19 sec on Intel P9500 (64bit Linux, 4 GByte RAM)
###########################################################
## hgu95a
system.time(eset.hgu95a.pmonly <- robloxbioc(spikein.hgu95a, bg.correct = TRUE, 
                                             pmcorrect = FALSE, normalize = TRUE, 
                                             add.constant = 0))
eset.hgu95a.log2.pmonly <- eset.hgu95a.pmonly
exprs(eset.hgu95a.log2.pmonly) <- log2(exprs(eset.hgu95a.pmonly))
roblox.hgu95a.pmonly <- assessSpikeIn(eset.hgu95a.log2.pmonly, method.name = "roblox (PM)")
roblox.hgu95a.pmonly.2 <- assessSpikeIn2(eset.hgu95a.log2.pmonly, method.name = "roblox (PM)")

## hgu133a
system.time(eset.hgu133a.pmonly <- robloxbioc(spikein.hgu133a, bg.correct = TRUE, 
                                             pmcorrect = FALSE, normalize = TRUE, 
                                             add.constant = 0))
eset.hgu133a.log2.pmonly <- eset.hgu133a.pmonly
exprs(eset.hgu133a.log2.pmonly) <- log2(exprs(eset.hgu133a.pmonly))
roblox.hgu133a.pmonly <- assessSpikeIn(eset.hgu133a.log2.pmonly, method.name = "roblox (PM)")
roblox.hgu133a.pmonly.2 <- assessSpikeIn2(eset.hgu133a.log2.pmonly, method.name = "roblox (PM)")


###########################################################
## Example 4: MAS 5.0 + 32 with PM only
## both computations take about 19 sec on Intel P9500 (64bit Linux, 4 GByte RAM)
###########################################################
## hgu95a
system.time(eset.hgu95a.pmonly32 <- robloxbioc(spikein.hgu95a, bg.correct = TRUE, 
                                             pmcorrect = FALSE, normalize = TRUE, 
                                             add.constant = 32))
eset.hgu95a.log2.pmonly32 <- eset.hgu95a.pmonly32
exprs(eset.hgu95a.log2.pmonly32) <- log2(exprs(eset.hgu95a.pmonly32))
roblox.hgu95a.pmonly32 <- assessSpikeIn(eset.hgu95a.log2.pmonly32, method.name = "roblox + 32 (PM)")
roblox.hgu95a.pmonly32.2 <- assessSpikeIn2(eset.hgu95a.log2.pmonly32, method.name = "roblox + 32 (PM)")

## hgu133a
system.time(eset.hgu133a.pmonly32 <- robloxbioc(spikein.hgu133a, bg.correct = TRUE, 
                                             pmcorrect = FALSE, normalize = TRUE, 
                                             add.constant = 32))
eset.hgu133a.log2.pmonly32 <- eset.hgu133a.pmonly32
exprs(eset.hgu133a.log2.pmonly32) <- log2(exprs(eset.hgu133a.pmonly32))
roblox.hgu133a.pmonly32 <- assessSpikeIn(eset.hgu133a.log2.pmonly32, method.name = "roblox + 32 (PM)")
roblox.hgu133a.pmonly32.2 <- assessSpikeIn2(eset.hgu133a.log2.pmonly32, method.name = "roblox + 32 (PM)")


###############################################################################
## Figure 1: The MA plot shows log fold change as a function of mean log 
## expression level. A set of 14 arrays representing a single experiment from 
## the Affymetrix spike-in data are used for this plot. A total of 13 sets of 
## fold changes are generated by comparing the first array in the set to each 
## of the others. Genes are symbolized by numbers representing the nominal 
## log2 fold change for the gene. Non-differentially expressed genes with 
## observed fold changes larger than 2 are plotted in red. All other probesets 
## are represented with black dots.
###############################################################################
## hgu95a
par(mfrow = c(3, 2))
affycompPlot(roblox.hgu95a$MA)
affycompPlot(roblox.hgu95a32$MA)
affycompPlot(roblox.hgu95a.pmonly$MA)
affycompPlot(roblox.hgu95a.pmonly32$MA)
affycompPlot(mas5.assessment$MA)
affycompPlot(rma.assessment$MA)

## hgu133a
par(mfrow = c(3, 2))
affycompPlot(roblox.hgu133a$MA)
affycompPlot(roblox.hgu133a32$MA)
affycompPlot(roblox.hgu133a.pmonly$MA)
affycompPlot(roblox.hgu133a.pmonly32$MA)
affycompPlot(mas5.assessment.133$MA)
affycompPlot(rma.assessment.133$MA)


###############################################################################
## Figure 4a: Average observed log2 intensity plotted against nominal log2 
## concentration for each spiked-in gene for all arrays in Affymetrix spike-In 
## experiment
###############################################################################
## hgu95a
par(mfrow = c(3, 2))
affycomp.figure4a(roblox.hgu95a$Signal)
affycomp.figure4a(roblox.hgu95a32$Signal)
affycomp.figure4a(roblox.hgu95a.pmonly$Signal)
affycomp.figure4a(roblox.hgu95a.pmonly32$Signal)
affycomp.figure4a(mas5.assessment$Signal)
affycomp.figure4a(rma.assessment$Signal)

## hgu133a
par(mfrow = c(3, 2))
affycomp.figure4a(roblox.hgu133a$Signal)
affycomp.figure4a(roblox.hgu133a32$Signal)
affycomp.figure4a(roblox.hgu133a.pmonly$Signal)
affycomp.figure4a(roblox.hgu133a.pmonly32$Signal)
affycomp.figure4a(mas5.assessment.133$Signal)
affycomp.figure4a(rma.assessment.133$Signal)

## Comparison plot
## hgu95a
affycomp.compfig4a(list(roblox.hgu95a$Signal, 
                        roblox.hgu95a32$Signal,
                        roblox.hgu95a.pmonly$Signal,
                        roblox.hgu95a.pmonly32$Signal,
                        mas5.assessment$Signal, 
                        rma.assessment$Signal), 
                   method.names = c("roblox", "roblox + 32", "roblox (PM)", 
                                    "roblox + 32 (PM)", "MAS 5.0", "RMA"))

## hgu133a
affycomp.compfig4a(list(roblox.hgu133a$Signal, 
                        roblox.hgu133a32$Signal,
                        roblox.hgu133a.pmonly$Signal,
                        roblox.hgu133a.pmonly32$Signal,
                        mas5.assessment.133$Signal, 
                        rma.assessment.133$Signal), 
                   method.names = c("roblox", "roblox + 32", "roblox (PM)", 
                                    "roblox + 32 (PM)", "MAS 5.0", "RMA"))


###############################################################################
## Figure 5: A typical identification rule for differential expression filters 
## genes with fold change exceeding a given threshold. This figure shows 
## average ROC curves which offer a graphical representation of both 
## specificity and sensitivity for such a detection rule. 
## a) Average ROC curves based on comparisons with nominal fold changes 
## ranging from 2 to 4096. 
## b) As a) but with nominal fold changes equal to 2.
###############################################################################
## Figure 5a:
## hgu95a
par(mfrow = c(3, 2))
affycomp.figure5a(roblox.hgu95a$FC)
affycomp.figure5a(roblox.hgu95a32$FC)
affycomp.figure5a(roblox.hgu95a.pmonly$FC)
affycomp.figure5a(roblox.hgu95a.pmonly32$FC)
affycomp.figure5a(mas5.assessment$FC)
affycomp.figure5a(rma.assessment$FC)

## hgu133a
par(mfrow = c(3, 2))
affycomp.figure5a(roblox.hgu133a$FC)
affycomp.figure5a(roblox.hgu133a32$FC)
affycomp.figure5a(roblox.hgu133a.pmonly$FC)
affycomp.figure5a(roblox.hgu133a.pmonly32$FC)
affycomp.figure5a(mas5.assessment.133$FC)
affycomp.figure5a(rma.assessment.133$FC)

## Comparison plot
## hgu95a
affycomp.compfig5a(list(roblox.hgu95a$FC,
                        mas5.assessment$FC), 
                   method.names = c("roblox", "MAS 5.0"))
affycomp.compfig5a(list(roblox.hgu95a$FC,
                        roblox.hgu95a32$FC,
                        roblox.hgu95a.pmonly$FC,
                        roblox.hgu95a.pmonly32$FC,
                        mas5.assessment$FC, 
                        rma.assessment$FC), 
                   method.names = c("roblox", "roblox + 32", "roblox (PM)", 
                                    "roblox + 32 (PM)", "MAS 5.0", "RMA"))

## hgu133a
affycomp.compfig5a(list(roblox.hgu133a$FC,
                        mas5.assessment.133$FC), 
                   method.names = c("roblox", "MAS 5.0"))
affycomp.compfig5a(list(roblox.hgu133a$FC, 
                        roblox.hgu133a32$FC,
                        roblox.hgu133a.pmonly$FC,
                        roblox.hgu133a.pmonly32$FC,
                        mas5.assessment.133$FC, 
                        rma.assessment.133$FC), 
                   method.names = c("roblox", "roblox + 32", "roblox (PM)", 
                                    "roblox + 32 (PM)", "MAS 5.0", "RMA"))

## Figure 5b:
## hgu95a
par(mfrow = c(3, 2))
affycomp.figure5b(roblox.hgu95a$FC)
affycomp.figure5b(roblox.hgu95a32$FC)
affycomp.figure5b(roblox.hgu95a.pmonly$FC)
affycomp.figure5b(roblox.hgu95a.pmonly32$FC)
affycomp.figure5b(mas5.assessment$FC)
affycomp.figure5b(rma.assessment$FC)

## hgu133a
par(mfrow = c(3, 2))
affycomp.figure5b(roblox.hgu133a$FC)
affycomp.figure5b(roblox.hgu133a32$FC)
affycomp.figure5b(roblox.hgu133a.pmonly$FC)
affycomp.figure5b(roblox.hgu133a.pmonly32$FC)
affycomp.figure5b(mas5.assessment.133$FC)
affycomp.figure5b(rma.assessment.133$FC)

## Comparison plot
## hgu95a
affycomp.compfig5b(list(roblox.hgu95a$FC,
                        mas5.assessment$FC), 
                   method.names = c("roblox", "MAS 5.0"))
affycomp.compfig5b(list(roblox.hgu95a$FC, 
                        roblox.hgu95a32$FC,
                        roblox.hgu95a.pmonly$FC,
                        roblox.hgu95a.pmonly32$FC,
                        mas5.assessment$FC, 
                        rma.assessment$FC), 
                   method.names = c("roblox", "roblox + 32", "roblox (PM)", 
                                    "roblox + 32 (PM)", "MAS 5.0", "RMA"))

## hgu133a
affycomp.compfig5b(list(roblox.hgu133a$FC,
                        mas5.assessment.133$FC), 
                   method.names = c("roblox", "MAS 5.0"))
affycomp.compfig5b(list(roblox.hgu133a$FC, 
                        roblox.hgu133a32$FC,
                        roblox.hgu133a.pmonly$FC,
                        roblox.hgu133a.pmonly32$FC,
                        mas5.assessment.133$FC, 
                        rma.assessment.133$FC), 
                   method.names = c("roblox", "roblox + 32", "roblox (PM)", 
                                    "roblox + 32 (PM)", "MAS 5.0", "RMA"))


###############################################################################
## Figure 6:            
## a) Observed log fold changes plotted against nominal log fold changes. The
## dashed lines represent highest, 25th highest, 100th highest, 25th 
## percentile, 75th percentile, smallest 100th, smallest 25th, and smallest 
## log fold change for the genes that were not differentially expressed. 
## b) Like a) but the observed fold changes were calculated for spiked in 
## genes with nominal concentrations no higher than 2pM.
###############################################################################
## Figure 6a:
## hgu95a
par(mfrow = c(3, 2))
affycomp.figure6a(roblox.hgu95a$FC)
affycomp.figure6a(roblox.hgu95a32$FC)
affycomp.figure6a(roblox.hgu95a.pmonly$FC)
affycomp.figure6a(roblox.hgu95a.pmonly32$FC)
affycomp.figure6a(mas5.assessment$FC)
affycomp.figure6a(rma.assessment$FC)

## hgu133a
par(mfrow = c(3, 2))
affycomp.figure6a(roblox.hgu133a$FC)
affycomp.figure6a(roblox.hgu133a32$FC)
affycomp.figure6a(roblox.hgu133a.pmonly$FC)
affycomp.figure6a(roblox.hgu133a.pmonly32$FC)
affycomp.figure6a(mas5.assessment.133$FC)
affycomp.figure6a(rma.assessment.133$FC)

## Figure 6b:
## hgu95a
par(mfrow = c(3, 2))
affycomp.figure6b(roblox.hgu95a$FC)
affycomp.figure6b(roblox.hgu95a32$FC)
affycomp.figure6b(roblox.hgu95a.pmonly$FC)
affycomp.figure6b(roblox.hgu95a.pmonly32$FC)
affycomp.figure6b(mas5.assessment$FC)
affycomp.figure6b(rma.assessment$FC)

## hgu133a
par(mfrow = c(3, 2))
affycomp.figure6b(roblox.hgu133a$FC)
affycomp.figure6b(roblox.hgu133a32$FC)
affycomp.figure6b(roblox.hgu133a.pmonly$FC)
affycomp.figure6b(roblox.hgu133a.pmonly32$FC)
affycomp.figure6b(mas5.assessment.133$FC)
affycomp.figure6b(rma.assessment.133$FC)


###############################################################################
## Table
###############################################################################
## hgu95a
tab.hgu95a <- tableAll(roblox.hgu95a, mas5.assessment)
round(tab.hgu95a, 4)
tab.hgu95a <- tableAll(roblox.hgu95a, roblox.hgu95a32, 
                        roblox.hgu95a.pmonly, roblox.hgu95a.pmonly32, 
                        mas5.assessment, rma.assessment)
round(tab.hgu95a, 4)


system.time(mas5.res <- mas5(spikein.hgu95a))
system.time(rma.res <- rma(spikein.hgu95a))
mas5.exprs <- exprs(mas5.res)
exprs(mas5.res) <- log2(mas5.exprs)
mas5.ass <- assessSpikeIn2(mas5.res, method.name = "MAS 5.0")
rma.ass <- assessSpikeIn2(rma.res, method.name = "RMA")

tab.hgu95a.2 <- tableAll(roblox.hgu95a.2, mas5.ass)
round(tab.hgu95a.2, 4)
tab.hgu95a.2 <- tableAll(roblox.hgu95a.2, roblox.hgu95a32.2, 
                         roblox.hgu95a.pmonly.2, roblox.hgu95a.pmonly32.2,
                         mas5.ass, rma.ass)
round(tab.hgu95a.2, 4)


## hgu133a
tab.hgu133a <- tableAll(roblox.hgu133a, mas5.assessment.133)
round(tab.hgu133a, 4)
tab.hgu133a <- tableAll(roblox.hgu133a, roblox.hgu133a32, 
                        roblox.hgu133a.pmonly, roblox.hgu133a.pmonly32, 
                        mas5.assessment.133, rma.assessment.133)
round(tab.hgu133a, 4)

system.time(mas5.res.133 <- mas5(spikein.hgu133a))
system.time(rma.res.133 <- rma(spikein.hgu133a))
mas5.exprs <- exprs(mas5.res.133)
exprs(mas5.res.133) <- log2(mas5.exprs)
mas5.ass.133 <- assessSpikeIn2(mas5.res.133, method.name = "MAS 5.0")
rma.ass.133 <- assessSpikeIn2(rma.res.133, method.name = "RMA")

tab.hgu133a.2 <- tableAll(roblox.hgu133a.2, mas5.ass.133)
round(tab.hgu133a.2, 4)
tab.hgu133a.2 <- tableAll(roblox.hgu133a.2, roblox.hgu133a32.2, 
                        roblox.hgu133a.pmonly.2, roblox.hgu133a.pmonly32.2,
                        mas5.ass.133, rma.ass.133)
round(tab.hgu133a.2, 4)


## smaller table, more informative ...
## affycompTable does not work due to missing dilution data
## hgu95a
tab.hgu95a.small <- tab.hgu95a[c(1,2,6:8,15:17,9:11), ]
tab.hgu95a.small <- cbind(tab.hgu95a.small, "whatsgood" = c(1, 1, 1, 0, 16, 1, 0, 16, 0, 1, 1))
round(tab.hgu95a.small, 4)

## hgu133a
tab.hgu133a.small <- tab.hgu133a[c(1,2,6:8,15:17,9:11), ]
tab.hgu133a.small <- cbind(tab.hgu133a.small, "whatsgood" = c(1, 1, 1, 0, NA, 1, 0, NA, 0, 1, 1))
round(tab.hgu133a.small, 4)


## Figure 3 in Kohl and Deigner (2010)
res.sd <- assessSpikeInSD(eset.hgu133a.log2)
res.sd.133 <- assessSpikeInSD(eset.hgu133a.log2)
res.sd.mas <- assessSpikeInSD(mas5.res)
res.sd.mas.133 <- assessSpikeInSD(mas5.res.133)
#setEPS(width = 9, height = 6)
#postscript(file = "Figure3.eps")
par(mfrow = c(1, 2))
plot(res.sd.mas$xsmooth, res.sd.mas$ysmooth, type = "l", xlab = "mean log expression",
     ylab = "mean SD", main = "HGU95A", lwd = 2, ylim = c(0, max(res.sd.mas$ysmooth)),
     xlim = c(0, max(res.sd.mas$xsmooth)),
     panel.first = abline(h = seq(0, 1.2, by = 0.2), v = c(0, 5, 10, 15), lty = 2, col = "grey"))
lines(res.sd$xsmooth, res.sd$ysmooth, type = "l", lty = 2, lwd = 2)
legend("topright", c("biweight", "rmx"), lty = 1:2, lwd = 2, bg = "white")
plot(res.sd.133$xsmooth, res.sd.133$ysmooth, type = "l", lty = 2, lwd = 2, 
     xlab = "mean log expression", ylab = "mean SD", main = "HGU133A", 
     ylim = c(0, max(res.sd.mas$ysmooth)), 
     xlim = c(0, max(res.sd.133$xsmooth)),
     panel.first = abline(h = seq(0, 1.2, by = 0.2), v = c(0, 5, 10, 15), lty = 2, col = "grey"))
lines(res.sd.mas.133$xsmooth, res.sd.mas.133$ysmooth, type = "l", lwd = 2)
legend("topright", c("biweight", "rmx"), lty = 1:2, lwd = 2, bg = "white")
#dev.off()


## Table 4 in Kohl and Deigner (2010)
round(tableAll(roblox.hgu95a.2, mas5.ass)[c(1:3, 8:11),], 3)
round(tableAll(roblox.hgu133a.2, mas5.ass.133)[c(1:3, 8:11),], 3)

