### R code from vignette source 'trackeR.Rnw'

###################################################
### code chunk number 1: preliminaries
###################################################
options(width = 70, prompt = "R> ", continue = "+  ")
library("trackeR")
library("ggplot2")
set.seed(403)
cache <- FALSE


###################################################
### code chunk number 2: trackeR.Rnw:1362-1367
###################################################
filepath <- system.file("extdata", "2013-06-08-090442.TCX", package = "trackeR")
## read raw data
runDF <- readTCX(file = filepath, timezone = "GMT")
## runDF is a data.frame with the following structure
str(runDF)


###################################################
### code chunk number 3: trackeR.Rnw:2014-2016
###################################################
## turn the raw data in runDF into a trackeRdata object
runTr0 <- trackeRdata(runDF)


###################################################
### code chunk number 4: trackeR.Rnw:2024-2026
###################################################
runTr1 <- readContainer(filepath, type = "tcx", timezone = "GMT")
identical(runTr0, runTr1)


###################################################
### code chunk number 5: dataLoad
###################################################
data(run, package = "trackeR")
data(runs, package = "trackeR")


###################################################
### code chunk number 6: trackeR.Rnw:2067-2068 (eval = FALSE)
###################################################
## plot(runs, session = 1:3)


###################################################
### code chunk number 7: defaultPlot
###################################################
plot(runs, session = 1:3)


###################################################
### code chunk number 8: trackeR.Rnw:2092-2093 (eval = FALSE)
###################################################
## plotRoute(runs, session = 4, zoom = 13, source = "osm")


###################################################
### code chunk number 9: routePlot
###################################################
plotRoute(runs, session = 4, zoom = 13, source = "osm")


###################################################
### code chunk number 10: summarySessionsPrint
###################################################
summary(runs, session = 1:2)
runSummary <- summary(runs, session = 1)
print(runSummary, digits = 3)


###################################################
### code chunk number 11: summarySessions (eval = FALSE)
###################################################
## runSummaryFull <- summary(runs)
## plot(runSummaryFull, group = c("total", "moving"),
##     what = c("avgSpeed", "distance", "duration", "avgHeartRate"))


###################################################
### code chunk number 12: summarySessionsPlot
###################################################
runSummaryFull <- summary(runs)
plot(runSummaryFull, group = c("total", "moving"),
    what = c("avgSpeed", "distance", "duration", "avgHeartRate"))


###################################################
### code chunk number 13: zones (eval = FALSE)
###################################################
## runZones <- zones(runs[1:4], what = "speed", breaks = list(speed = c(0, 2:6, 12.5)))
## ## if breaks is a named list, argument 'what' can be left unspecified
## runZones <- zones(runs[1:4], breaks = list(speed = c(0, 2:6, 12.5)))
## ## if only a single variable is to be evaluated, 'breaks' can also be a vector
## runZones <- zones(runs[1:4], what = "speed", breaks = c(0, 2:6, 12.5))
## plot(runZones)


###################################################
### code chunk number 14: zonesPlot
###################################################
runZones <- zones(runs[1:4], what = "speed", breaks = list(speed = c(0, 2:6, 12.5)))
## if breaks is a named list, argument 'what' can be left unspecified
runZones <- zones(runs[1:4], breaks = list(speed = c(0, 2:6, 12.5)))
## if only a single variable is to be evaluated, 'breaks' can also be a vector
runZones <- zones(runs[1:4], what = "speed", breaks = c(0, 2:6, 12.5))
plot(runZones)


###################################################
### code chunk number 15: Wprime (eval = FALSE)
###################################################
## wexp <- Wprime(runs, session = 11, quantity = "expended",
##                cp = 4, version = "2012")
## plot(wexp, scaled = TRUE)


###################################################
### code chunk number 16: WprimePlot2
###################################################
wexp <- Wprime(runs, session = 11, quantity = "expended", cp = 4, version = "2012")
plot(wexp, scaled = TRUE)


###################################################
### code chunk number 17: distrProfiles (eval = FALSE)
###################################################
## dProfile <- distributionProfile(runs, session = 1:4,
##     what = c("speed", "heart.rate"),
##     grid = list(speed = seq(0, 12.5, by = 0.05), heart.rate = seq(0, 250)))
## plot(dProfile, multiple = TRUE)


###################################################
### code chunk number 18: dProfile0 (eval = FALSE)
###################################################
## set.seed(1)
## dProfile <- distributionProfile(runs, session = 1:4,
##     what = c("speed", "heart.rate"),
##     grid = list(speed = seq(0, 12.5, by = 0.05), heart.rate = seq(0, 250)))
## dProfileS <- smoother(dProfile, cores = 2)


###################################################
### code chunk number 19: dProfile
###################################################
if(cache & file.exists("example_dProfile.rda")) {
  load("example_dProfile.rda")
} else {
set.seed(1)
dProfile <- distributionProfile(runs, session = 1:4,
    what = c("speed", "heart.rate"),
    grid = list(speed = seq(0, 12.5, by = 0.05), heart.rate = seq(0, 250)))
dProfileS <- smoother(dProfile, cores = 2)
if(cache) {
  save(dProfile, dProfileS, file = "example_dProfile.rda")
} else {
  if(file.exists("example_dProfile.rda")) file.remove("example_dProfile.rda")
}
}


###################################################
### code chunk number 20: dprofilePlot
###################################################
plot(dProfileS, smooth = FALSE, multiple = TRUE)


###################################################
### code chunk number 21: concentrationProfiles (eval = FALSE)
###################################################
## cProfile <- concentrationProfile(dProfile, what = "speed")
## plot(cProfile, multiple = TRUE)


###################################################
### code chunk number 22: cprofilePlot
###################################################
cProfile <- concentrationProfile(dProfile, what = "speed")
plot(cProfile, multiple = TRUE, cores = 2)
## cProfileS <- concentrationProfile(dProfileS)
## plot(cProfileS, what = "speed", smooth = FALSE, multiple = TRUE)


###################################################
### code chunk number 23: trackeR.Rnw:2597-2603
###################################################
## get the units for the variables in run
getUnits(run)

## change the unit of speed into miles per hour
runTr2 <- changeUnits(run, variable = "speed", unit = "mi_per_h")
getUnits(runTr2)


###################################################
### code chunk number 24: summarySessionsUnits
###################################################
m_per_s2ft_per_h <- function(x) x * 3937/1200 * 3600
changeUnits(runSummary, variable = "speed", unit = "ft_per_h")


###################################################
### code chunk number 25: thresholdPlots (eval = FALSE)
###################################################
## ## without thresholds
## plot(runs, session = 4, what = "speed", threshold = FALSE)
## ## with default thresholds
## plot(runs, session = 4, what = "speed")
## ## with default thresholds and smoothing
## plot(runs, session = 4, what = "speed", smooth = TRUE, fun = "median", width = 20)
## ## thresholding and smoothing outside of plot method
## run4 <- threshold(runs[4])
## run4S <- smoother(run4, what = "speed", fun = "median", width = 20)
## plot(run4S, what = "speed", smooth = FALSE) 


###################################################
### code chunk number 26: thresholdPlots1
###################################################
plot(runs, session = 4, what = "speed", threshold = FALSE)


###################################################
### code chunk number 27: thresholdPlots2
###################################################
plot(runs, session = 4, what = "speed") + ggplot2::expand_limits(y = c(0, 21))


###################################################
### code chunk number 28: thresholdPlots3
###################################################
plot(runs, session = 4, what = "speed", smooth = TRUE, fun = "median", width = 20, cores = 2) + ggplot2::expand_limits(y = c(0, 12.5))


###################################################
### code chunk number 29: thresholdPlots4
###################################################
run4 <- threshold(runs[4])
run4S <- smoother(run4, what = "speed", fun = "median", width = 20, cores = 2)
plot(run4S, what = "speed", smooth = FALSE) + ggplot2::expand_limits(y = c(0, 12.5))


###################################################
### code chunk number 30: AppData
###################################################
library("trackeR")
data(runs, package = "trackeR")
runsSummary <- summary(runs)


###################################################
### code chunk number 31: AppProfiles0 (eval = FALSE)
###################################################
## library("trackeR")
## ## load data
## data(runs, package = "trackeR")
## ## apply default thresholds
## runsT <- threshold(runs)
## ## get and smooth distribution profiles 
## dpRuns <- distributionProfile(runsT, what = "speed")
## dpRunsS <- smoother(dpRuns, cores = 2)
## ## get concentration profiles
## cpRuns <- concentrationProfile(dpRunsS)


###################################################
### code chunk number 32: AppProfiles
###################################################
if(cache & file.exists("appProfiles.rda")) {
    load("appProfiles.rda")
} else {
library("trackeR")
## load data
data(runs, package = "trackeR")
## apply default thresholds
runsT <- threshold(runs)
## get and smooth distribution profiles 
dpRuns <- distributionProfile(runsT, what = "speed")
dpRunsS <- smoother(dpRuns, cores = 2)
## get concentration profiles
cpRuns <- concentrationProfile(dpRunsS)
if(cache) {
   save(dpRuns, dpRunsS, cpRuns, file = "appProfiles.rda")
} else {
    if(file.exists("appProfiles.rda")) file.remove("appProfiles.rda")
}
}


###################################################
### code chunk number 33: AppCPplot
###################################################
plot(cpRuns, multiple = TRUE, smooth = FALSE) + theme(legend.position="none")


###################################################
### code chunk number 34: AppFunPrep
###################################################
## prepare functional data
library("fda")
gridSpeed <- seq(0, 12.5, length = 251)
sp <- matrix(unlist(cpRuns$speed), ncol = 250, byrow = TRUE,
             dimnames = list(names(cpRuns$speed), gridSpeed[-1]))
spfd <- Data2fd(argvals = gridSpeed[-1], y = t(sp))
## fit functional PCA
sppca <- pca.fd(spfd, nharm = 4)
## share of variance
varprop <- round(sppca$varprop * 100); names(varprop) <- 1:4
varprop


###################################################
### code chunk number 35: mypcaplot
###################################################
mypcaplot <- function (x, nx = 128, pointplot = TRUE, harm = 0, expand = 0, 
    cycle = FALSE, xlab = "argvals", ...) 
{
    pcafd <- x
    if (!(inherits(pcafd, "pca.fd"))) 
        stop("Argument 'x' is not a pca.fd object.")
    harmfd <- pcafd[[1]]
    basisfd <- harmfd$basis
    rangex <- basisfd$rangeval
    if (length(nx) > 1) {
        argvals <- nx
        nx <- length(x)
    }
    else {
        argvals <- seq(rangex[1], rangex[2], length = nx)
    }
    fdmat <- eval.fd(argvals, harmfd)
    meanmat <- eval.fd(argvals, pcafd$meanfd)
    dimfd <- dim(fdmat)
    nharm <- dimfd[2]
    plotsPerPg <- sum(par("mfrow"))
    harm <- as.vector(harm)
    if (harm[1] == 0) 
        harm <- (1:nharm)
    if (length(dimfd) == 2) {
        for (jharm in 1:length(harm)) {
            if (jharm == 2) {
                op <- par(ask = TRUE)
                on.exit(par(op))
            }
            iharm <- harm[jharm]
            if (expand == 0) {
                fac <- sqrt(pcafd$values[iharm])
            }
            else {
                fac <- expand
            }
            vecharm <- fdmat[, iharm]
            pcmat <- cbind(meanmat + fac * vecharm, meanmat - 
                fac * vecharm)
            if (pointplot) 
                plottype <- "p"
            else plottype <- "l"
            percentvar <- round(100 * pcafd$varprop[iharm], 1)
            plot(argvals, meanmat, type = "l", ylim = c(min(pcmat), 
                max(pcmat)), ylab = paste("Harmonic", iharm), xlab = xlab,
                main = paste("PCA Function", iharm, "(Percentage of Variability", 
                  percentvar, ")"))
            if (pointplot) {
                points(argvals, pcmat[, 1], pch = "+")
                points(argvals, pcmat[, 2], pch = "-")
            }
            else {
                lines(argvals, pcmat[, 1], lty = 2)
                lines(argvals, pcmat[, 2], lty = 3)
            }
        }
    }
    else {
        if (cycle && dimfd[3] == 2) {
            meanmat <- drop(meanmat)
            for (jharm in 1:length(harm)) {
                if (jharm == 2) {
                  op <- par(ask = TRUE)
                  on.exit(par(op))
                }
                iharm <- harm[jharm]
                {
                  if (expand == 0) 
                    fac <- 2 * sqrt(pcafd$values[iharm])
                  else fac <- expand
                }
                matharm <- fdmat[, iharm, ]
                mat1 <- meanmat + fac * matharm
                mat2 <- meanmat - fac * matharm
                if (pointplot) 
                  plottype <- "p"
                else plottype <- "l"
                percentvar <- round(100 * pcafd$varprop[iharm], 
                  1)
                plot(meanmat[, 1], meanmat[, 2], type = plottype, 
                  xlim = c(min(c(mat1[, 1], mat2[, 1])), max(c(mat1[, 
                    1], mat2[, 1]))), ylim = c(min(c(mat1[, 2], 
                    mat2[, 2])), max(c(mat1[, 2], mat2[, 2]))), 
                  main = paste("PCA Function", iharm, "(Percentage of Variability", 
                    percentvar, ")"), ...)
                if (pointplot) {
                  points(mat1[, 1], mat1[, 2], pch = "+")
                  points(mat2[, 1], mat2[, 2], pch = "-")
                }
                else {
                  lines(mat1[, 1], mat1[, 2], lty = 2)
                  lines(mat2[, 1], mat2[, 2], lty = 3)
                }
            }
        }
        else {
            for (jharm in 1:length(harm)) {
                if (jharm == 2) {
                  op <- par(ask = TRUE)
                  on.exit(par(op))
                }
                iharm <- harm[jharm]
                fac <- {
                  if (expand == 0) 
                    sqrt(pcafd$values[iharm])
                  else expand
                }
                meanmat <- drop(meanmat)
                matharm <- fdmat[, iharm, ]
                nvar <- dim(matharm)[2]
                for (jvar in 1:nvar) {
                  pcmat <- cbind(meanmat[, jvar] + fac * matharm[, 
                    jvar], meanmat[, jvar] - fac * matharm[, 
                    jvar])
                  if (pointplot) 
                    plottype <- "p"
                  else plottype <- "l"
                  percentvar <- round(100 * pcafd$varprop[iharm], 
                    1)
                  plot(argvals, meanmat[, jvar], type = plottype, xlab = xlab,
                    ylab = paste("Harmonic", iharm), sub = paste("PCA Function", 
                      iharm, "(Percentage of Variability", percentvar, 
                      ")"), main = dimnames(fdmat)[[3]][jvar], 
                    ...)
                  if (pointplot) {
                    points(argvals, pcmat[, 1], pch = "+")
                    points(argvals, pcmat[, 2], pch = "-")
                  }
                  else {
                    lines(argvals, pcmat[, 1], lty = 2)
                    lines(argvals, pcmat[, 2], lty = 3)
                  }
                }
            }
        }
    }
    invisible(NULL)
}


###################################################
### code chunk number 36: AppHarmonicsPlotSpeed
###################################################
##source("mypcaplot.R")
pp <- FALSE ## dotted = -, dashed = +
par(mfrow = c(2,1), ask = FALSE)
mypcaplot(sppca, harm = 1:2, pointplot = pp, xlab = "Speed [m/s]")
par(mfrow = c(1,1), ask = TRUE)


###################################################
### code chunk number 37: AppScoresPlot1
###################################################
## plot scores vs summary statistics
scoresSP <- data.frame(sppca$scores)
names(scoresSP) <- paste0("speed_pc", 1:4)
d <- cbind(runsSummary, scoresSP)

library("ggplot2")
## pc1 ~ session duration (moving)
d$durationMoving <- as.numeric(d$durationMoving)
ggplot(d) + geom_point(aes(x = durationMoving, y = speed_pc1)) + theme_bw() + labs(x = "duration moving [min]", y = "PC1")


###################################################
### code chunk number 38: AppScoresPlot2
###################################################
## pc2 ~ avg speed (moving)
ggplot(d) + geom_point(aes(x = avgSpeedMoving, y = speed_pc2)) + theme_bw() + labs(x = "average speed moving [m/s]", y = "PC2")


