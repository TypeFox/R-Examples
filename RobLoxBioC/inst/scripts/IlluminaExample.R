###############################################################################
## Illumina Example
###############################################################################

###############################################################################
## References:
## Dunning, M.J., Smith, M.L., Ritchie, M.E., Tavare, S.:
## beadarray: R classes and methods for Illumina bead-based data. 
## Bioinformatics 2007, 23(16):2183-4.
##
## M.J. Dunning, N.L. Barbosa-Morais, A.G. Lynch, S. Tavaré and M.E. Ritchie:
## Statistical issues in the analysis of Illumina data.
## BMC Bioinformatics 2008, 9:85.
###############################################################################

###############################################################################
## Data:
## Can be obtained via
## http://www.compbio.group.cam.ac.uk/Resources/spike/index.html
###############################################################################

## Load the required packages
library(beadarray)
library(gplots)
library(RColorBrewer)
library(RobLoxBioC)

###############################################################################
## Extract all *.zip file to directory "SpikeInData".
## Copy spike_targets.txt to directory "SpikeInData".
##
## Code to read the bead level data from the directory "SpikeInData"
##
## NOTE: reading in the raw data for the entire experiment requires at
## least 4Gb of RAM for each processing method.  
###############################################################################

###########################################################
## Read targets information
targets <- read.table("./SpikeInData/spike_targets.txt",header=TRUE)
arraynms <- as.character(targets$ArrayNo)

## Use sharpened, subtracted data from text files
spikeInData <- readIllumina(path = "./SpikeInData", arrayNames=arraynms[1:2], 
                            useImages=FALSE, textType=".csv")
#save(spikeInData, compress = TRUE, file = "spikeInData.RData")
#load(file = "spikeInData.RData")

## takes about 9 hours on Intel P9500 (64bit Linux, 4 GByte RAM)
system.time(minKD.Illumina <- KolmogorovMinDist(spikeInData, Norm(), imagesPerArray = 2))
save(minKD.Illumina, compress = TRUE, file = "minKD_Illumina.RData")

## takes about 9 hours on Intel P9500 (64bit Linux, 4 GByte RAM)
system.time(minKD.Illumina.log <- KolmogorovMinDist(spikeInData, Norm(), log = TRUE, imagesPerArray = 2))
save(minKD.Illumina.log, compress = TRUE, file = "minKD_Illumina_log.RData")

## load the results from R-forge ...
con <- url("http://robast.r-forge.r-project.org/data/minKD_Illumina.RData")
load(file = con)
close(con)
con <- url("http://robast.r-forge.r-project.org/data/minKD_Illumina_log.RData")
load(file = con)
close(con)

## takes more than 90 min on Intel P9500 (64bit Linux, 4 GByte RAM)
ns <- c(10:70)
M <- length(ns)
minKD.Illumina.norm <- matrix(NA, nrow = 50000, ncol = M)
colnames(minKD.Illumina.norm) <- ns
#for(i in seq_len(M)){
for(i in 27:61){
    tm <- proc.time()
    print(ns[i])
    temp <- matrix(rnorm(50000*ns[i]), ncol = ns[i])
    minKD.Illumina.norm[,i] <- KolmogorovMinDist(temp, Norm())$dist
    cat("Dauer:\t", proc.time()-tm, "\n")
    save(minKD.Illumina.norm, compress = TRUE, file = "minKD_Illumina_norm1.RData")
}

## load the results from R-forge
con <- url("http://robast.r-forge.r-project.org/data/minKD_Illumina_norm.RData")
load(file = con)
close(con)

#######################################
## Figure 5 in Kohl and Deigner (2010)
#######################################
res1 <- split(as.vector(minKD.Illumina$dist), as.vector(minKD.Illumina$n))[30:50]
res2 <- split(as.vector(minKD.Illumina.log$dist), as.vector(minKD.Illumina.log$n))[30:50]
res3 <- lapply(as.data.frame(minKD.Illumina.norm[,21:41]), function(x) x)
uni.n <- rep(30:50, 3)

#setEPS(height = 6, width = 9)
#postscript(file = "Figure5.eps")
par(mar = c(4, 4, 3, 1))
plot(0, 0, type = "n", ylim = c(-0.01, 0.4), xlim = c(0.5, 65.5), 
     panel.first = abline(h = seq(0, 0.35, by = 0.05), lty = 2, col = "grey"), 
     main = "Minimum Kolmogorov distance", 
     ylab = "minimum Kolmogorov distance", 
     xlab = "sample size", axes = FALSE)
axis(1, c(1:21, 23:43, 45:65), labels = uni.n, cex.axis = 0.6)
axis(2, seq(0, 0.35, by = 0.05), labels = seq(0, 0.35, by = 0.05), las = 2,
     cex.axis = 0.8)
box()
boxplot(c(res1, res2, res3), at = c(1:21, 23:43, 45:65), add = TRUE, pch = 20, 
        names = FALSE, axes = FALSE)
abline(h = c(0.055, 0.093), lty = 2, lwd = 1.5)
abline(v = c(22, 44), lwd = 1.5)
text(c(10, 33, 55), rep(0.38, 3), labels = c("Bead Level Data", "log Bead Level Data", "Normal Samples"),
     font = 2)
lines(1:21, 1/(2*(30:50)), lwd = 2)
lines(23:43, 1/(2*(30:50)), lwd = 2)
lines(45:65, 1/(2*(30:50)), lwd = 2)
legend("bottomleft", legend = "minimal possible distance", lty = 1, 
       bg = "white", cex = 0.8)
#dev.off()

## Comparison of quantiles
## Figure 6 in Kohl and Deigner (2010)
res1 <- split(as.vector(minKD.Illumina$dist), as.vector(minKD.Illumina$n))[15:65]
res2 <- split(as.vector(minKD.Illumina.log$dist), as.vector(minKD.Illumina.log$n))[15:65]
res3 <- lapply(as.data.frame(minKD.Illumina.norm), function(x) x)[6:56]

#setEPS(height = 6, width = 9)
#postscript(file = "Figure6.eps")
par(mar = c(4, 4, 3, 1))
plot(15:65, sapply(res3, quantile, prob = 0.99), type = "l", lwd = 2, xlab = "sample size", 
     ylab = "quantile of mimimum Kolmogorov distances",
     main = "50% and 99% quantiles of minimum Kolmogorov distances", 
     ylim = c(0.05, 0.23),
     panel.first = abline(h = c(0.05, 0.1, 0.15, 0.2), v = seq(15, 65, by = 5), 
                          lty = 2, col = "grey"))
lines(15:65, sapply(res1, quantile, prob = 0.99), lwd = 2, lty = 2)
lines(15:65, sapply(res2, quantile, prob = 0.99), lwd = 2, lty = 3)
lines(15:65, sapply(res3, quantile, prob = 0.5), lwd = 2, lty = 1)
lines(15:65, sapply(res1, quantile, prob = 0.5), lwd = 2, lty = 2)
lines(15:65, sapply(res2, quantile, prob = 0.5), lwd = 2, lty = 3)
text(22, 0.18, "99% quantiles", font = 2)
text(22, 0.115, "50% quantiles", font = 2)
legend("topright", legend = c("normal samples", "bead level data", "log bead level data"),
       lty = 1:3, lwd = 2, bg = "white")
#dev.off()


#load(file = "spikeInData.RData")
## takes about 100 sec on Intel P9500 (64bit Linux, 4 GByte RAM)
system.time(res.ill <- createBeadSummaryData(spikeInData, log = TRUE, imagesPerArray = 2))
## takes about 500 sec on Intel P9500 (64bit Linux, 4 GByte RAM)
system.time(res.rmx <- robloxbioc(spikeInData, imagesPerArray = 2))

## cf. function assessSpikeInSD of package affycomp
library(splines)
assessSpikeInSD <- function(BSdata, genenames, method.name = NULL, span = 1/3){
    spikein <- match(genenames, featureNames(BSdata))
    y <- apply(exprs(BSdata)[-spikein, ], 1, sd, na.rm = TRUE)
    x <- rowMeans(exprs(BSdata)[-spikein, ], na.rm = TRUE)
    smooth1 <- loess(y ~ x, span = 1/3, family = "gaussian", degree = 1)
    x2 <- sort(x)[seq(1, length(x), length = 100)]
    y2 <- smooth1$fitted[order(x)][seq(1, length(x), length = 100)]
    list(x = x, y = y, xsmooth = x2, ysmooth = y2, loess = smooth1, 
         method.name = method.name, what = "SpikeInSD")
}

spikecsv <- read.csv("spikeins_profile.csv")
genenames <- spikecsv[,"ProbeID"]
ill.SD <- assessSpikeInSD(res.ill, genenames = genenames, method.name = "Illumina")
rmx.SD <- assessSpikeInSD(res.rmx, genenames = genenames, method.name = "rmx estimator")

## Figure 7 in Kohl and Deigner (2010)
setEPS(height = 6, width = 9)
postscript(file = "Figure7.eps")
plot(ill.SD$xsmooth, ill.SD$ysmooth, type = "l", xlab = "mean log expression",
     ylab = "mean SD", main = "Spike-in data of Dunning et al. (2008)", lwd = 2,
     panel.first = abline(h = c(0.1, 0.12, 0.14, 0.16), v = seq(6, 16, 2), lty = 2, col = "grey"))
lines(rmx.SD$xsmooth, rmx.SD$ysmooth, type = "l", lty = 2, lwd = 2)
legend("topright", c("Illumina", "rmx"), lty = 1:2, lwd = 2, bg = "white")
dev.off()


## Table 6 in Kohl and Deigner (2009)
quantile(ill.SD$y, prob = c(0.25, 0.5, 0.75, 0.99))
quantile(rmx.SD$y, prob = c(0.25, 0.5, 0.75, 0.99))

## cf. functions assessMA2 and tableMA2 of package affycomp
genenames <- spikecsv[,"ProbeID"]
assessMA2 <- function(BSdata, genenames){
    NCOMP <- choose(12, 2)*4
    I <- c(1:6, 25:30, 7:12, 31:36, 13:18, 37:42, 19:24, 43:48)
    spikein <- match(genenames, featureNames(BSdata))
    mat <- exprs(BSdata)[,I]
    quants <- matrix(0, nrow(mat) - length(spikein), NCOMP)
    count <- 0
    for (i in c(1:11, 13:23, 25:35, 37:47)) {
        for (j in (i + 1):(12 * ceiling(i/12))) {
            count <- count + 1
            fc <- mat[, j] - mat[, i]
            quants[, count] <- sort(fc[-spikein], na.last = TRUE)
        }
    }
    nulls <- rowMedians(quants)
    temp <- quantile(nulls, prob = c(0.99, 0.999), na.rm = TRUE)
    c("null log-fc IQR" = IQR(nulls, na.rm = TRUE), 
      "null log-fc 99%" = temp[1],
      "null log-fc 99.9%" = temp[2])
}
assessMA2(BSdata = res.ill, genenames = genenames)
assessMA2(BSdata = res.rmx, genenames = genenames)


## other procedures
system.time(res.ill1 <- createBeadSummaryData(spikeInData, log = TRUE, imagesPerArray = 2,
                                              method = "mean"))
system.time(res.ill2 <- createBeadSummaryData(spikeInData, log = TRUE, imagesPerArray = 2,
                                              method = "median"))
system.time(res.ill3 <- createBeadSummaryData(spikeInData, log = TRUE, imagesPerArray = 2,
                                              method = "trim"))
system.time(res.ill4 <- createBeadSummaryData(spikeInData, log = TRUE, imagesPerArray = 2,
                                              method = "winsorize"))
system.time(res.rmx1 <- robloxbioc(spikeInData, imagesPerArray = 2, eps.upper = 0.025))
system.time(res.rmx2 <- robloxbioc(spikeInData, imagesPerArray = 2, eps.lower = 0.01))
ill1.SD <- assessSpikeInSD(res.ill1, genenames = genenames, method.name = "mean")
ill2.SD <- assessSpikeInSD(res.ill2, genenames = genenames, method.name = "median")
ill3.SD <- assessSpikeInSD(res.ill3, genenames = genenames, method.name = "trim")
ill4.SD <- assessSpikeInSD(res.ill4, genenames = genenames, method.name = "winsorize")
rmx1.SD <- assessSpikeInSD(res.rmx1, genenames = genenames, method.name = "rmx estimator")
rmx2.SD <- assessSpikeInSD(res.rmx2, genenames = genenames, method.name = "rmx estimator")
quantile(ill.SD$y, prob = c(0.25, 0.5, 0.75, 0.99))
quantile(ill1.SD$y, prob = c(0.25, 0.5, 0.75, 0.99))
quantile(ill2.SD$y, prob = c(0.25, 0.5, 0.75, 0.99))
quantile(ill3.SD$y, prob = c(0.25, 0.5, 0.75, 0.99))
quantile(ill4.SD$y, prob = c(0.25, 0.5, 0.75, 0.99))
quantile(rmx.SD$y, prob = c(0.25, 0.5, 0.75, 0.99))
quantile(rmx1.SD$y, prob = c(0.25, 0.5, 0.75, 0.99))
quantile(rmx2.SD$y, prob = c(0.25, 0.5, 0.75, 0.99))
assessMA2(BSdata = res.ill, genenames = genenames)
assessMA2(BSdata = res.ill1, genenames = genenames)
assessMA2(BSdata = res.ill2, genenames = genenames)
assessMA2(BSdata = res.ill3, genenames = genenames)
assessMA2(BSdata = res.ill4, genenames = genenames)
assessMA2(BSdata = res.rmx, genenames = genenames)
assessMA2(BSdata = res.rmx1, genenames = genenames)
assessMA2(BSdata = res.rmx2, genenames = genenames)


###############################################################################
## The following example is based on the R code of Mark Dunning and Matt Ritchie
## available under
## http://www.compbio.group.cam.ac.uk/Resources/spike/scripts/Analysis.R
##
## This file was slightly adapted and code for the computation of 
## rmx-estimators was added.
###############################################################################

########################################
## Analysis of Illumina in Spike data set
##
## November 2007
##
## Mark Dunning and Matt Ritchie
########################################

###############################################################################
## Extract all *.zip file to directory "SpikeInData".
## Copy spike_targets.txt to directory "SpikeInData".
##
## Code to read the bead level data from the directory "SpikeInData"
##
## NOTE: reading in the raw data for the entire experiment requires at
## least 4Gb of RAM for each processing method.  
## For this reason the data is read in sequentially saved, then removed
## before moving on to the next method 
###############################################################################


###########################################################
## Modified Figure 2 of Dunning et al (2008) 
## plot results from simulation study (varying number of outliers 
## were introduced at the saturation level to assess how many outliers can be 
## tolerated before serious bias is introduced in to the ## expression measures)

targets <- read.table("./SpikeInData/spike_targets.txt", header=TRUE, sep=" ")
arraynms <- as.character(targets$ArrayNo)
narrays <- length(arraynms)

## introduce outliers in 2nd BeadChip, at varying %
## Use sharpened, subtracted data from text files
bld.sharpen.array2 <- readIllumina(path = "./SpikeInData", arrayNames=arraynms[13:24], 
                                   useImages=FALSE, textType=".csv")

## percentage of beads to be saturated
per.out <- c(0, 0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4)

nchips <- 12
arraynms <- arrayNames(bld.sharpen.array2)

## Saturated
ill.summ.array2 <- list(NoBeads=list(), exprs=list(), se.exprs=list()) 
ill.summ2.array2 <- med.summ.array2 <- mean.summ.array2 <- trim.summ.array2 <- win.summ.array2 <- ill.summ.array2
rmx.summ.array2 <- ill.summ.array2

set.seed(06092007)
for(l in 1:length(per.out)) {
    cat(l, "of 9\n")
    simdataarray2 <- copyBeadLevelList(bld.sharpen.array2)
    if(per.out[l]!=0) {
        for(i in 1:nchips) {
            cat(i, "of 12\n")
        ## no gross error model!
#            nsamp <- round(per.out[l]*numBeads(simdataarray2, array=i),0)
#            ind <- sample(1:numBeads(simdataarray2, array=i), nsamp)
        ## gross error model!
#            ind <- as.logical(rbinom(numBeads(simdataarray2, array=i), prob = per.out[l], size = 1))
#            ind <- (1:numBeads(simdataarray2, array=i))[ind]
#            simdataarray2@beadData[[arraynms[i]]]$G[ind] <- 2^16
        ## problem: there are bead types where more than 50% of the values are contaminated
        ## => there is no meaningful estimator which can handle this!
        ## Hence, we contiminate bead-type wise
            sel <- which(simdataarray2@beadData[[arraynms[i]]]$ProbeID != 0)
            pr <- simdataarray2@beadData[[arraynms[i]]]$ProbeID[sel]
            probes <- sort(unique(pr))
            indices <- NULL
            for(j in seq(along = probes)){
                ind <- pr == probes[j]
                ## < 50% of values contaminated
                repeat{
                    cont <- as.logical(rbinom(sum(ind), prob = per.out[l], size = 1))
                    if(sum(cont) < sum(ind)/2) break
                }
                indices <- c(indices, sel[ind][cont])
            }
            simdataarray2@beadData[[arraynms[i]]]$G[indices] <- 2^16
        }
    }
    tmp.ill <- createBeadSummaryData(simdataarray2, method="illumina", log=FALSE, n=3, imagesPerArray=2)
    tmp.ill2 <- createBeadSummaryData(simdataarray2, method="illumina", log=FALSE, n=2, imagesPerArray=2) 
    tmp.med <- createBeadSummaryData(simdataarray2, method="median", log=FALSE, imagesPerArray=2)
    tmp.mean <- createBeadSummaryData(simdataarray2, method="mean", log=FALSE, imagesPerArray=2)
    tmp.trim <- createBeadSummaryData(simdataarray2, method="trim", trim=0.1, log=FALSE, imagesPerArray=2)
    tmp.win <- createBeadSummaryData(simdataarray2, method="winsorize", trim=0.1, log=FALSE, imagesPerArray=2)
    tmp.rmx <- robloxbioc(simdataarray2, imagesPerArray=2, eps = max(per.out[l], 0.05))
    ill.summ.array2$NoBeads[[l]] = NoBeads(tmp.ill)
    ill.summ.array2$exprs[[l]] = exprs(tmp.ill)
    ill.summ.array2$se.exprs[[l]] = se.exprs(tmp.ill)
    ill.summ2.array2$NoBeads[[l]] = NoBeads(tmp.ill2)
    ill.summ2.array2$exprs[[l]] = exprs(tmp.ill2)
    ill.summ2.array2$se.exprs[[l]] = se.exprs(tmp.ill2)
    med.summ.array2$NoBeads[[l]] = NoBeads(tmp.med)
    med.summ.array2$exprs[[l]] = exprs(tmp.med)
    med.summ.array2$se.exprs[[l]] = se.exprs(tmp.med)
    mean.summ.array2$NoBeads[[l]] = NoBeads(tmp.mean)
    mean.summ.array2$exprs[[l]] = exprs(tmp.mean)
    mean.summ.array2$se.exprs[[l]] = se.exprs(tmp.mean)
    win.summ.array2$NoBeads[[l]] = NoBeads(tmp.win)
    win.summ.array2$exprs[[l]] = exprs(tmp.win)
    win.summ.array2$se.exprs[[l]] = se.exprs(tmp.win)
    trim.summ.array2$NoBeads[[l]] = NoBeads(tmp.trim)
    trim.summ.array2$exprs[[l]] = exprs(tmp.trim)
    trim.summ.array2$se.exprs[[l]] = se.exprs(tmp.trim)
    rmx.summ.array2$NoBeads[[l]] = NoBeads(tmp.rmx)
    rmx.summ.array2$exprs[[l]] = exprs(tmp.rmx)
    rmx.summ.array2$se.exprs[[l]] = se.exprs(tmp.rmx)
    rm(simdataarray2, tmp.ill, tmp.med, tmp.mean, tmp.trim, tmp.win, tmp.rmx)
}
save(ill.summ.array2, ill.summ2.array2, med.summ.array2, mean.summ.array2, 
     win.summ.array2, trim.summ.array2, rmx.summ.array2, file="sim.summary.bias.raw.rda")

load(file="sim.summary.bias.raw.rda")

## calculate bias - assume original, complete data set gave true values
bias.ill.array2 <- list(NoBeads=list(), exprs=list(), se.exprs=list())
bias.ill2.array2 <- bias.med.array2 <- bias.mean.array2 <- bias.trim.array2 <- bias.win.array2 <- bias.ill.array2
bias.rmx.array2 <- bias.ill.array2
for(l in 1:length(per.out)) {
    cat(l, "\n")
    bias.ill.array2$exprs[[l]] <- (ill.summ.array2$exprs[[l]]-mean.summ.array2$exprs[[1]])^2
    bias.ill.array2$se.exprs[[l]] <- ill.summ.array2$se.exprs[[l]]^2*ill.summ.array2$NoBeads[[l]]
    bias.ill.array2$NoBeads[[l]] <- mean.summ.array2$NoBeads[[1]]-ill.summ.array2$NoBeads[[l]]

    bias.ill2.array2$exprs[[l]] <- (ill.summ2.array2$exprs[[l]]-mean.summ.array2$exprs[[1]])^2
    bias.ill2.array2$se.exprs[[l]] <- ill.summ2.array2$se.exprs[[l]]^2*ill.summ2.array2$NoBeads[[l]]
    bias.ill2.array2$NoBeads[[l]] <- mean.summ.array2$NoBeads[[1]]-ill.summ2.array2$NoBeads[[l]]

    bias.med.array2$exprs[[l]] <- (med.summ.array2$exprs[[l]]-mean.summ.array2$exprs[[1]])^2
    bias.med.array2$se.exprs[[l]] <- med.summ.array2$se.exprs[[l]]^2*med.summ.array2$NoBeads[[l]]
    bias.med.array2$NoBeads[[l]] <- mean.summ.array2$NoBeads[[1]]-med.summ.array2$NoBeads[[l]]
  
    bias.mean.array2$exprs[[l]] <- (mean.summ.array2$exprs[[l]]-mean.summ.array2$exprs[[1]])^2
    bias.mean.array2$se.exprs[[l]] <- mean.summ.array2$se.exprs[[l]]^2*mean.summ.array2$NoBeads[[l]]
    bias.mean.array2$NoBeads[[l]] <- mean.summ.array2$NoBeads[[1]]-mean.summ.array2$NoBeads[[l]]

    bias.trim.array2$exprs[[l]] <- (trim.summ.array2$exprs[[l]]-mean.summ.array2$exprs[[1]])^2
    bias.trim.array2$se.exprs[[l]] <- trim.summ.array2$se.exprs[[l]]^2*trim.summ.array2$NoBeads[[l]]
    bias.trim.array2$NoBeads[[l]] <- mean.summ.array2$NoBeads[[1]]-trim.summ.array2$NoBeads[[l]]

    bias.win.array2$exprs[[l]] <- (win.summ.array2$exprs[[l]]-mean.summ.array2$exprs[[1]])^2
    bias.win.array2$se.exprs[[l]] <- win.summ.array2$se.exprs[[l]]^2*win.summ.array2$NoBeads[[l]]
    bias.win.array2$NoBeads[[l]] <- mean.summ.array2$NoBeads[[1]]-win.summ.array2$NoBeads[[l]]

    bias.rmx.array2$exprs[[l]] <- (rmx.summ.array2$exprs[[l]]-mean.summ.array2$exprs[[1]])^2
    bias.rmx.array2$se.exprs[[l]] <- rmx.summ.array2$se.exprs[[l]]^2*rmx.summ.array2$NoBeads[[l]]
    bias.rmx.array2$NoBeads[[l]] <- mean.summ.array2$NoBeads[[1]]-rmx.summ.array2$NoBeads[[l]]
}

## calculate average bias from `replicate arrays
avebias.ill.array2 <- avebias.ill2.array2 <- avebias.mean.array2 <- avebias.med.array2 <- avebias.trim.array2 <- avebias.win.array2 <- NULL
avebias.rmx.array2 <- NULL
for(l in 1:length(per.out)) {
    avebias.ill.array2[l] <- mean(as.vector(bias.ill.array2$exprs[[l]]), na.rm=TRUE)
    avebias.ill2.array2[l] <- mean(as.vector(bias.ill2.array2$exprs[[l]]), na.rm=TRUE)
    avebias.med.array2[l] <- mean(as.vector(bias.med.array2$exprs[[l]]), na.rm=TRUE)
    avebias.mean.array2[l] <- mean(as.vector(bias.mean.array2$exprs[[l]]), na.rm=TRUE)
    avebias.win.array2[l] <- mean(as.vector(bias.win.array2$exprs[[l]]), na.rm=TRUE)
    avebias.trim.array2[l] <- mean(as.vector(bias.trim.array2$exprs[[l]]), na.rm=TRUE)
    avebias.rmx.array2[l] <- mean(as.vector(bias.rmx.array2$exprs[[l]]), na.rm=TRUE)
}

## calculate variance
avevar.ill.array2 <- avevar.ill2.array2 <- avevar.mean.array2 <- avevar.med.array2 <- avevar.trim.array2 <- avevar.win.array2 <- NULL
avevar.rmx.array2 <- NULL
for(l in 1:length(per.out)) {
    avevar.ill.array2[l] <- mean(as.vector(bias.ill.array2$se.exprs[[l]]), na.rm=TRUE)
    avevar.ill2.array2[l] <- mean(as.vector(bias.ill2.array2$se.exprs[[l]]), na.rm=TRUE)
    avevar.med.array2[l] <- mean(as.vector(bias.med.array2$se.exprs[[l]]), na.rm=TRUE)
    avevar.mean.array2[l] <- mean(as.vector(bias.mean.array2$se.exprs[[l]]), na.rm=TRUE)
    avevar.win.array2[l] <- mean(as.vector(bias.win.array2$se.exprs[[l]]), na.rm=TRUE)
    avevar.trim.array2[l] <- mean(as.vector(bias.trim.array2$se.exprs[[l]]), na.rm=TRUE)
    avevar.rmx.array2[l] <- mean(as.vector(bias.rmx.array2$se.exprs[[l]]), na.rm=TRUE)
}

## calculate MSE
avemse.ill.array2 <- avemse.ill2.array2 <- avemse.mean.array2 <- avemse.med.array2 <- avemse.trim.array2 <- avemse.win.array2 <- NULL
avemse.rmx.array2 <- NULL
for(l in 1:length(per.out)) {
    avemse.ill.array2[l] <- mean(as.vector(bias.ill.array2$se.exprs[[l]])+as.vector(bias.ill.array2$exprs[[l]]), na.rm=TRUE)
    avemse.ill2.array2[l] <- mean(as.vector(bias.ill2.array2$se.exprs[[l]])+as.vector(bias.ill2.array2$exprs[[l]]), na.rm=TRUE)
    avemse.med.array2[l] <- mean(as.vector(bias.med.array2$se.exprs[[l]])+as.vector(bias.med.array2$exprs[[l]]), na.rm=TRUE)
    avemse.mean.array2[l] <- mean(as.vector(bias.mean.array2$se.exprs[[l]])+as.vector(bias.mean.array2$exprs[[l]]), na.rm=TRUE)
    avemse.win.array2[l] <- mean(as.vector(bias.win.array2$se.exprs[[l]])+as.vector(bias.win.array2$exprs[[l]]), na.rm=TRUE)
    avemse.trim.array2[l] <- mean(as.vector(bias.trim.array2$se.exprs[[l]])+as.vector(bias.trim.array2$exprs[[l]]), na.rm=TRUE)
    avemse.rmx.array2[l] <- mean(as.vector(bias.rmx.array2$se.exprs[[l]])+as.vector(bias.rmx.array2$exprs[[l]]), na.rm=TRUE)
}

lwd <- 2
sel <- 1:9
x <- per.out[sel]*100
myCol <- brewer.pal(5, "Set1")

pdf("Figure2.pdf", width=8, height=5)
par(mfrow=c(1,3))
par(mar=c(3,4,1.5,0.15), oma=c(1,0,0,0))
plot(x, log2(avebias.ill.array2[sel]), lwd=lwd, type="l", xlab="", ylab=expression(log[2]("average bias^2")), 
     main="Bias^2", ylim=c(0,30), col = myCol[1], 
     panel.first = abline(v = c(0, 10, 20, 30, 40), h = seq(0, 30, by = 5), lty = 2, col = "grey"))
points(x, log2(avebias.med.array2[sel]), type="l", col=myCol[2], lwd=lwd)
#points(x, avebias.trim.array2[sel], type="l", col=myCol[3], lwd=lwd)
points(x, log2(avebias.mean.array2[sel]), type="l", col=myCol[4], lwd=lwd)
points(x, log2(avebias.rmx.array2[sel]), type="l", col=myCol[5], lwd=lwd)
legend("bottomright",legend=c("Illumina", "median", "mean", "rmx"), col=myCol[c(1,2,4,5)], lwd=2)
#plot(x, numout.ill.array2[sel], lwd=lwd, type="l", xlab="", ylab="% observations removed by summary method", main="(b)", ylim=c(0,40))
#points(x, numout.med.array2[sel], type="l", col=2, lwd=lwd)
#points(x, numout.trim.array2[sel], type="l", col=3, lwd=lwd)
#points(x, numout.mean.array2[sel], type="l", col=4, lwd=lwd)
#abline(0,1,col="gray", lty=2)
plot(x, log2(avevar.ill.array2)[sel], lwd=lwd, type="l", xlab="", ylab=expression(log[2]("average variance")), 
     main="Variance", ylim=c(0,30),col = myCol[1], 
     panel.first = abline(v = c(0, 10, 20, 30, 40), h = seq(0, 30, by = 5), lty = 2, col = "grey"))
points(x, log2(avevar.med.array2)[sel], type="l", col=myCol[2], lwd=lwd)
#points(x, log2(avevar.trim.array2)[sel], type="l", col=myCol[3], lwd=lwd)
points(x, log2(avevar.mean.array2)[sel], type="l", col=myCol[4], lwd=lwd)
points(x, log2(avevar.rmx.array2)[sel], type="l", col=myCol[5], lwd=lwd)
plot(x, log2(avemse.ill.array2)[sel], lwd=lwd, type="l", xlab="", ylab=expression(log[2]("average MSE")), 
     main="MSE", ylim=c(0,30),col = myCol[1], 
     panel.first = abline(v = c(0, 10, 20, 30, 40), h = seq(0, 30, by = 5), lty = 2, col = "grey"))
points(x, log2(avemse.med.array2)[sel], type="l", col=myCol[2], lwd=lwd)
#points(x, log2(avemse.trim.array2)[sel], type="l", col=myCol[3], lwd=lwd)
points(x, log2(avemse.mean.array2)[sel], type="l", col=myCol[4], lwd=lwd)
points(x, log2(avemse.rmx.array2)[sel], type="l", col=myCol[5], lwd=lwd)
mtext("% outliers simulated", side=1, outer=TRUE)
dev.off()


## log-scale
## Saturated
ill.log.summ.array2 <- list(NoBeads=list(), exprs=list(), se.exprs=list()) 
ill.log.summ2.array2 <- med.log.summ.array2 <- mean.log.summ.array2 <- trim.log.summ.array2 <- win.log.summ.array2 <- ill.log.summ.array2
ill.log.rmx.array2 <- ill.log.summ.array2

set.seed(06092007)
for(l in 1:length(per.out)) {
    cat(l, "\n")
    simdataarray2 <- copyBeadLevelList(bld.sharpen.array2)
    if(per.out[l]!=0) {
        for(i in 1:nchips) {
            cat(i, "of 12\n")
        ## no gross error model!
#            nsamp <- round(per.out[l]*numBeads(simdataarray2, array=i),0)
#            ind <- sample(1:numBeads(simdataarray2, array=i), nsamp)
        ## gross error model!
#            ind <- as.logical(rbinom(numBeads(simdataarray2, array=i), prob = per.out[l], size = 1))
#            ind <- (1:numBeads(simdataarray2, array=i))[ind]
#            simdataarray2@beadData[[arraynms[i]]]$G[ind] <- 2^16
        ## problem: there are bead types where more than 50% of the values are contaminated
        ## => there is no meaningful estimator which can handle this!
            sel <- which(simdataarray2@beadData[[arraynms[i]]]$ProbeID != 0)
            pr <- simdataarray2@beadData[[arraynms[i]]]$ProbeID[sel]
            probes <- sort(unique(pr))
            indices <- NULL
            for(j in seq(along = probes)){
                ind <- pr == probes[j]
                ## < 50% of values contaminated
                repeat{
                    cont <- as.logical(rbinom(sum(ind), prob = per.out[l], size = 1))
                    if(sum(cont) < sum(ind)/2) break
                }
                indices <- c(indices, sel[ind][cont])
            }
            simdataarray2@beadData[[arraynms[i]]]$G[indices] <- 2^16
        }
    }
    tmp.ill.log <- createBeadSummaryData(simdataarray2, method="illumina", log=TRUE, n=3, imagesPerArray=2)
    tmp.ill2.log <- createBeadSummaryData(simdataarray2, method="illumina", log=TRUE, n=2, imagesPerArray=2) 
    tmp.med.log <- createBeadSummaryData(simdataarray2, method="median", log=TRUE, imagesPerArray=2)
    tmp.mean.log <- createBeadSummaryData(simdataarray2, method="mean", log=TRUE, imagesPerArray=2)
    tmp.trim.log <- createBeadSummaryData(simdataarray2, method="trim", trim=0.1, log=TRUE, imagesPerArray=2)
    tmp.win.log <- createBeadSummaryData(simdataarray2, method="winsorize", trim=0.1, log=TRUE, imagesPerArray=2)
    tmp.rmx.log <- robloxbioc(simdataarray2, log = TRUE, imagesPerArray=2, eps = max(per.out[l], 0.05))
    ill.log.summ.array2$NoBeads[[l]] <- NoBeads(tmp.ill.log)
    ill.log.summ.array2$exprs[[l]] <- exprs(tmp.ill.log)
    ill.log.summ.array2$se.exprs[[l]] <- se.exprs(tmp.ill.log)
    ill.log.summ2.array2$NoBeads[[l]] <- NoBeads(tmp.ill2.log)
    ill.log.summ2.array2$exprs[[l]] <- exprs(tmp.ill2.log)
    ill.log.summ2.array2$se.exprs[[l]] <- se.exprs(tmp.ill2.log)
    med.log.summ.array2$NoBeads[[l]] <- NoBeads(tmp.med.log)
    med.log.summ.array2$exprs[[l]] <- exprs(tmp.med.log)
    med.log.summ.array2$se.exprs[[l]] <- se.exprs(tmp.med.log)
    mean.log.summ.array2$NoBeads[[l]] <- NoBeads(tmp.mean.log)
    mean.log.summ.array2$exprs[[l]] <- exprs(tmp.mean.log)
    mean.log.summ.array2$se.exprs[[l]] <- se.exprs(tmp.mean.log)
    win.log.summ.array2$NoBeads[[l]] <- NoBeads(tmp.win.log)
    win.log.summ.array2$exprs[[l]] <- exprs(tmp.win.log)
    win.log.summ.array2$se.exprs[[l]] <- se.exprs(tmp.win.log)
    trim.log.summ.array2$NoBeads[[l]] <- NoBeads(tmp.trim.log)
    trim.log.summ.array2$exprs[[l]] <- exprs(tmp.trim.log)
    trim.log.summ.array2$se.exprs[[l]] <- se.exprs(tmp.trim.log)
    rmx.log.summ.array2$NoBeads[[l]] <- NoBeads(tmp.rmx.log)
    rmx.log.summ.array2$exprs[[l]] <- exprs(tmp.rmx.log)
    rmx.log.summ.array2$se.exprs[[l]] <- se.exprs(tmp.rmx.log)
    rm(simdataarray2, tmp.ill.log, tmp.med.log, tmp.mean.log, tmp.trim.log, tmp.win.log, tmp.rmx.log)
}
save(ill.log.summ.array2, ill.log.summ2.array2, med.log.summ.array2, mean.log.summ.array2, 
     win.log.summ.array2, trim.log.summ.array2, rmx.log.summ.array2, file="sim.summary.bias.log.rda")

bias.ill.log.array2 <- list(NoBeads=list(), exprs=list(), se.exprs=list())
bias.ill2.log.array2 <- bias.med.log.array2 <- bias.mean.log.array2 <- bias.trim.log.array2 <- bias.win.log.array2 <- bias.ill.log.array2
bias.rmx.log.array2 <- bias.ill.log.array2
for(l in 1:length(per.out)) {
    cat(l, "\n")
    bias.ill.log.array2$exprs[[l]] <- ill.log.summ.array2$exprs[[l]]-mean.log.summ.array2$exprs[[1]]
    bias.ill.log.array2$se.exprs[[l]] <- ill.log.summ.array2$se.exprs[[l]]^2*ill.log.summ.array2$NoBeads[[l]]
    bias.ill.log.array2$NoBeads[[l]] <- mean.log.summ.array2$NoBeads[[1]]-ill.log.summ.array2$NoBeads[[l]]

    bias.ill2.log.array2$exprs[[l]] <- ill.log.summ2.array2$exprs[[l]]-mean.log.summ.array2$exprs[[1]]
    bias.ill2.log.array2$se.exprs[[l]] <- ill.log.summ2.array2$se.exprs[[l]]^2*ill.log.summ2.array2$NoBeads[[l]]
    bias.ill2.log.array2$NoBeads[[l]] <- mean.log.summ.array2$NoBeads[[1]]-ill.log.summ2.array2$NoBeads[[l]]

    bias.med.log.array2$exprs[[l]] <- med.log.summ.array2$exprs[[l]]-mean.log.summ.array2$exprs[[1]]
    bias.med.log.array2$NoBeads[[l]] <- mean.log.summ.array2$NoBeads[[1]]-med.log.summ.array2$NoBeads[[l]]
    
    bias.mean.log.array2$exprs[[l]] <- mean.log.summ.array2$exprs[[l]]-mean.log.summ.array2$exprs[[1]]
    bias.mean.log.array2$se.exprs[[l]] <- mean.log.summ.array2$se.exprs[[l]]^2*mean.log.summ.array2$NoBeads[[l]]
    bias.mean.log.array2$NoBeads[[l]] <- mean.log.summ.array2$NoBeads[[1]]-mean.log.summ.array2$NoBeads[[l]]

    bias.trim.log.array2$exprs[[l]] <- trim.log.summ.array2$exprs[[l]]-mean.log.summ.array2$exprs[[1]]
    bias.trim.log.array2$se.exprs[[l]] <- trim.log.summ.array2$se.exprs[[l]]^2*trim.log.summ.array2$NoBeads[[l]]
    bias.trim.log.array2$NoBeads[[l]] <- mean.log.summ.array2$NoBeads[[1]]-trim.log.summ.array2$NoBeads[[l]]

    bias.win.log.array2$exprs[[l]] <- win.log.summ.array2$exprs[[l]]-mean.log.summ.array2$exprs[[1]]
    bias.win.log.array2$se.exprs[[l]] <- win.log.summ.array2$se.exprs[[l]]^2*win.log.summ.array2$NoBeads[[l]]
    bias.win.log.array2$NoBeads[[l]] <- mean.log.summ.array2$NoBeads[[1]]-win.log.summ.array2$NoBeads[[l]]

    bias.rmx.log.array2$exprs[[l]] <- rmx.log.summ.array2$exprs[[l]]-mean.log.summ.array2$exprs[[1]]
    bias.rmx.log.array2$se.exprs[[l]] <- rmx.log.summ.array2$se.exprs[[l]]^2*rmx.log.summ.array2$NoBeads[[l]]
    bias.rmx.log.array2$NoBeads[[l]] <- mean.log.summ.array2$NoBeads[[1]]-rmx.log.summ.array2$NoBeads[[l]]
}

avebias.ill.log.array2 <- avebias.ill2.log.array2 <- avebias.mean.log.array2 <- avebias.med.log.array2 <- avebias.trim.log.array2 <- avebias.win.log.array2 <- NULL
avebias.rmx.log.array2 <- NULL
for(l in 1:length(per.out)) {
    avebias.ill.log.array2[l] <- mean(as.vector(bias.ill.log.array2$exprs[[l]]), na.rm=TRUE)
    avebias.ill2.log.array2[l] <- mean(as.vector(bias.ill2.log.array2$exprs[[l]]), na.rm=TRUE)
    avebias.med.log.array2[l] <- mean(as.vector(bias.med.log.array2$exprs[[l]]), na.rm=TRUE)
    avebias.mean.log.array2[l] <- mean(as.vector(bias.mean.log.array2$exprs[[l]]), na.rm=TRUE)
    avebias.win.log.array2[l] <- mean(as.vector(bias.win.log.array2$exprs[[l]]), na.rm=TRUE)
    avebias.trim.log.array2[l] <- mean(as.vector(bias.trim.log.array2$exprs[[l]]), na.rm=TRUE)
    avebias.rmx.log.array2[l] <- mean(as.vector(bias.rmx.log.array2$exprs[[l]]), na.rm=TRUE)
}

numout.ill.log.array2 <- numout.ill2.log.array2 <- numout.mean.log.array2 <- numout.med.log.array2 <- numout.trim.log.array2 <- numout.win.log.array2 <- NULL
numout.rmx.log.array2 <- NULL
for(l in 1:length(per.out)) {
    numout.ill.log.array2[l] <- mean(as.vector(bias.ill.log.array2$NoBeads[[l]]/mean.log.summ.array2$NoBeads[[1]]), na.rm=TRUE)*100
    numout.ill2.log.array2[l] <- mean(as.vector(bias.ill2.log.array2$NoBeads[[l]]/mean.log.summ.array2$NoBeads[[1]]), na.rm=TRUE)*100
    numout.med.log.array2[l] <- mean(as.vector(bias.med.log.array2$NoBeads[[l]]/mean.log.summ.array2$NoBeads[[1]]), na.rm=TRUE)*100
    numout.mean.log.array2[l] <- mean(as.vector(bias.mean.log.array2$NoBeads[[l]]/mean.log.summ.array2$NoBeads[[1]]), na.rm=TRUE)*100
    numout.win.log.array2[l] <- mean(as.vector(bias.win.log.array2$NoBeads[[l]]/mean.log.summ.array2$NoBeads[[1]]), na.rm=TRUE)*100
    numout.trim.log.array2[l] <- mean(as.vector(bias.trim.log.array2$NoBeads[[l]]/mean.log.summ.array2$NoBeads[[1]]), na.rm=TRUE)*100
    numout.rmx.log.array2[l] <- mean(as.vector(bias.rmx.log.array2$NoBeads[[l]]/mean.log.summ.array2$NoBeads[[1]]), na.rm=TRUE)*100
}

avevar.ill.log.array2 <- avevar.ill2.log.array2 <- avevar.mean.log.array2 <- avevar.med.log.array2 <- avevar.trim.log.array2 <- avevar.win.log.array2 <- NULL
avevar.rmx.log.array2 <- NULL
for(l in 1:length(per.out)) {
    avevar.ill.log.array2[l] <- mean(as.vector(bias.ill.log.array2$se.exprs[[l]]), na.rm=TRUE)
    avevar.ill2.log.array2[l] <- mean(as.vector(bias.ill2.log.array2$se.exprs[[l]]), na.rm=TRUE)
#    avevar.med.log.array2[l] <- mean(as.vector(bias.med.log.array2$se.exprs[[l]]), na.rm=TRUE)
    avevar.mean.log.array2[l] <- mean(as.vector(bias.mean.log.array2$se.exprs[[l]]), na.rm=TRUE)
    avevar.win.log.array2[l] <- mean(as.vector(bias.win.log.array2$se.exprs[[l]]), na.rm=TRUE)
    avevar.trim.log.array2[l] <- mean(as.vector(bias.trim.log.array2$se.exprs[[l]]), na.rm=TRUE)
    avevar.rmx.log.array2[l] <- mean(as.vector(bias.rmx.log.array2$se.exprs[[l]]), na.rm=TRUE)
}

pdf("Figure2b.pdf", width=12, height=8)
par(mfrow=c(1,2))
par(mar=c(3,4,1.5,0.15), oma=c(1,0,0,0))
plot(x, avebias.ill.log.array2[sel], lwd=lwd, type="l", xlab="", ylab="average bias", main="A", ylim=c(0,4))
points(x, avebias.med.log.array2[sel], type="l", col=2, lwd=lwd)
points(x, avebias.trim.log.array2[sel], type="l", col=3, lwd=lwd)
points(x, avebias.mean.log.array2[sel], type="l", col=4, lwd=lwd)
points(x, avebias.rmx.log.array2[sel], type="l", col=5, lwd=lwd)
legend("topleft",legend=c("Illumina", "median", "trimmed mean", "mean", "rmx"), col=1:5, lwd=2)
#plot(x, numout.ill.log.array2[sel], lwd=lwd, type="l", xlab="", ylab="% observations removed by summary method", main="(b)", ylim=c(0,40))
#points(x, numout.med.log.array2[sel], type="l", col=2, lwd=lwd)
#points(x, numout.trim.log.array2[sel], type="l", col=3, lwd=lwd)
#points(x, numout.mean.log.array2[sel], type="l", col=4, lwd=lwd)
#abline(0,1,col="gray", lty=2)
plot(x, log2(avevar.ill.log.array2)[sel], lwd=lwd, type="l", xlab="", ylab=expression(log[2]("average variance")), main="B", ylim=c(-5,5))
points(x, log2(avevar.trim.log.array2)[sel], type="l", col=3, lwd=lwd)
points(x, log2(avevar.mean.log.array2)[sel], type="l", col=4, lwd=lwd)
points(x, log2(avevar.rmx.log.array2)[sel], type="l", col=5, lwd=lwd)
mtext("% outliers simulated", side=1, outer=TRUE)
dev.off()
