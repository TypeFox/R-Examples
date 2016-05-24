#e###############################################################################
# # powerExample.R
# # R Versions: 2.15.0
# #
# # Author(s): Victor H. Cervantes
# #
# # Description: Shows how to use the package DFIT to calculate power
# #
# # Outputs: Example runs for power calculation in the DFIT approach for NCDIF
# #
# # File history:
# #   20120403: Creation
################################################################################

################################################################################
# # Load packages
################################################################################
library(DFIT)
#source("../R/IPR.R")
#source("../R/IRTSE.R")
#source("../R/MantelHaenszel.R")
#source("../R/RajuAreas.R")
#source("../R/NCDIF.R")
#library(ggplot2)

################################################################################
# # Load functions
################################################################################

################################################################################
# # Definition of input and output paths
################################################################################
#inPath  <- "../input/"
#outPath <- "../output/"
#srcPath <- "../src/"
#logPath <- "../log/"

################################################################################
# # Other global variables
################################################################################
nReplicates <- 3000
nFocal      <- 800
nReference  <- 2500
kRatio      <- nReference / nFocal

focalParam     <- list(mean = 0, sd = 1)
referenceParam <- list(mean = 0, sd = 1)


################################################################################
# # Uniform power example
################################################################################
itemParameters <- list(focal = cbind(rep(1, 51), seq(-0.5, 0.5, length = 51)),
                       reference = cbind(rep(1, 51), rep(0, 51)))

nullFocal      <- which(itemParameters[['focal']][, 2] ==
                        itemParameters[['reference']][, 2])

itemParametersNull <- lapply(itemParameters, function (x)
                             matrix(x[nullFocal, ], nrow = length(nullFocal)))
names(itemParametersNull) <- c('focal', 'reference')

# # True NCDIF
twoPlUniNcdifTrue    <- Ncdif(itemParameters, irtModel = '2pl',
                              focalDistribution = 'norm',
                              focalDistrExtra = focalParam,
                              logistic = FALSE) 


# Calculate cutoff point proposed algorithm
twoPlUniAse <- list()
twoPlUniAse[["focal"]] <- AseIrt(itemParameters = itemParameters[['focal']],
                                 distribution = 'norm',
                                 distributionParameters = focalParam,
                                 sampleSize = nFocal, irtModel = '2pl',
                                 logistic = FALSE)

twoPlUniAse[["reference"]] <- AseIrt(itemParameters =
                                     itemParameters[["reference"]],
                                     distribution = 'norm',
                                     distributionParameters = referenceParam,
                                     sampleSize = nReference, irtModel = "2pl",
                                     logistic = FALSE)

set.seed(29834328)


twoPlUniIpr          <- Ipr(itemParameters = itemParameters,
                            itemCovariances = twoPlUniAse,
                            nReplicates = nReplicates)

twoPlUniNcdif        <- IprNcdif(itemParameterList = twoPlUniIpr,
                                 irtModel = "2pl", logistic = FALSE,
                                 subdivisions = 1000)

cutoffPointEachSZUni <- CutoffIpr(quantiles = 0.95,
                                  iprStatistics =
                                  matrix(twoPlUniNcdif[nullFocal, ], nrow =
                                         length(nullFocal)))

# Calculate cutoff point current algorithm
set.seed(29834328)

twoPlUniAseCurrent                <- twoPlUniAse
twoPlUniAseCurrent[['focal']]     <- twoPlUniAseCurrent[['focal']][nullFocal]
twoPlUniAseCurrent[['reference']] <- lapply(twoPlUniAseCurrent[['reference']][nullFocal],
                                            '*', kRatio)

twoPlUniIprCurrent <- Ipr(itemParameters = itemParametersNull,
                          itemCovariances = twoPlUniAseCurrent,
                          nReplicates = nReplicates)

twoPlUniNcdifCurrent <- IprNcdif(itemParameterList = twoPlUniIprCurrent,
                                 irtModel = "2pl", logistic = FALSE,
                                 subdivisions = 1000)
cutoffPointUni <- CutoffIpr(quantiles = 0.95,
                            iprStatistics = matrix(twoPlUniNcdifCurrent,
                                                   nrow = length(nullFocal)))

# Calculate power
powerUni <- data.frame(Type = "Cutoff with only focal group",
                       b = itemParameters[['focal']][, 2],
                       NCDIF = (twoPlUniNcdifTrue *
                                sign(itemParameters[['focal']][, 2] - 
                                     itemParametersNull[['focal']][, 2])),
                       Power = rowMeans(twoPlUniNcdif >
                                        cutoffPointUni$quantiles))

powerUniEach <- data.frame(Type = "Cutoff with both groups",
                           b = itemParameters[['focal']][, 2],
                           NCDIF = (twoPlUniNcdifTrue *
                                    sign(itemParameters[['focal']][, 2] -
                                     itemParametersNull[['focal']][, 2])),
                           Power = rowMeans(twoPlUniNcdif >
                                            cutoffPointEachSZUni$quantiles))

power80P <- which(abs(powerUniEach[, 'Power'] - 0.8) ==
                  min(abs(powerUniEach[, 'Power'] - 0.8)))

powerUniJ <- rbind(powerUni, powerUniEach)

powerPlotUniDif <- ggplot(powerUniJ, aes(x = b, y = Power))
powerPlotUniDif <- powerPlotUniDif + xlim(-0.5, 0.5)
powerPlotUniDif <- powerPlotUniDif + geom_line(aes(colour = Type,
                                                   linetype = Type), size = 1.0)

powerPlotUniNcdif <- ggplot(powerUniJ, aes(x = NCDIF, y = Power))
powerPlotUniNcdif <- powerPlotUniNcdif + xlim(-0.012, 0.012)
powerPlotUniNcdif <- powerPlotUniNcdif + geom_line(aes(colour = Type,
                                                       linetype = Type),
                                                   size = 1.0)

powerPlotUniDif <- powerPlotUniDif + geom_hline(yintercept = 0.05,
                                                linetype = 'longdash',
                                                colour = 'grey')
powerPlotUniDif <- powerPlotUniDif + geom_hline(yintercept =
                                                powerUniEach[power80P, 'Power'],
                                                linetype = 'longdash',
                                                colour = 'grey')
powerPlotUniDif <- powerPlotUniDif + geom_hline(yintercept =
                                                powerUni[power80P, 'Power'],
                                                linetype = 'longdash',
                                                colour = 'grey')
powerPlotUniDif <- powerPlotUniDif + geom_vline(xintercept =
                                                powerUni[power80P, 'b'],
                                                linetype = 'longdash',
                                                colour = 'grey')
powerPlotUniDif <- powerPlotUniDif + geom_vline(xintercept =
                                                -powerUni[power80P, 'b'],
                                                linetype = 'longdash',
                                                colour = 'grey')
powerPlotUniDif


powerPlotUniNcdif <- powerPlotUniNcdif + geom_hline(yintercept = 0.05,
                                                    linetype = 'longdash',
                                                    colour = 'grey')
powerPlotUniNcdif <- powerPlotUniNcdif + geom_hline(yintercept =
                                                    powerUniEach[power80P, 'Power'],
                                                    linetype = 'longdash',
                                                    colour = 'grey')
powerPlotUniNcdif <- powerPlotUniNcdif + geom_hline(yintercept =
                                                    powerUni[power80P, 'Power'],
                                                    linetype = 'longdash',
                                                    colour = 'grey')
powerPlotUniNcdif <- powerPlotUniNcdif + geom_vline(xintercept =
                                                    powerUniEach[power80P, 'NCDIF'],
                                                    linetype = 'longdash',
                                                    colour = 'grey')
powerPlotUniNcdif <- powerPlotUniNcdif + geom_vline(xintercept =
                                                    -powerUniEach[power80P, 'NCDIF'],
                                                    linetype = 'longdash',
                                                    colour = 'grey')
powerPlotUniNcdif


################################################################################
# # NonUniform power example
################################################################################
itemParameters <- list(focal = cbind(unique(c(seq(0.5, 1, length = 26), 1,
                                              seq(1, 1 / 0.5, length = 26))), rep(0, 51)),
                       reference = cbind(rep(1, 51), rep(0, 51)))

nullFocal      <- which(itemParameters[['focal']][, 1] ==
                        itemParameters[['reference']][, 1])

itemParametersNull <- lapply(itemParameters, function (x)
                             matrix(x[nullFocal, ], nrow = length(nullFocal)))
names(itemParametersNull) <- c('focal', 'reference')

# # True NCDIF
twoPlNonUniNcdifTrue    <- Ncdif(itemParameters, irtModel = '2pl',
                              focalDistribution = 'norm',
                              focalDistrExtra = focalParam,
                              logistic = FALSE) 

# Calculate cutoff point proposed algorithm
twoPlNonUniAse <- list()
twoPlNonUniAse[["focal"]] <- AseIrt(itemParameters = itemParameters[['focal']],
                                    distribution = 'norm',
                                    distributionParameters = focalParam,
                                    sampleSize = nFocal, irtModel = '2pl',
                                    logistic = FALSE)

twoPlNonUniAse[["reference"]] <- AseIrt(itemParameters =
                                        itemParameters[["reference"]],
                                        distribution = 'norm',
                                        distributionParameters = referenceParam,
                                        sampleSize = nReference, irtModel = "2pl",
                                        logistic = FALSE)

set.seed(29834328)

twoPlNonUniIpr          <- Ipr(itemParameters = itemParameters,
                               itemCovariances = twoPlNonUniAse,
                               nReplicates = nReplicates)

twoPlNonUniNcdif        <- IprNcdif(itemParameterList = twoPlNonUniIpr,
                                    irtModel = "2pl", logistic = FALSE,
                                    subdivisions = 1000)

cutoffPointEachSZNonUni <- CutoffIpr(quantiles = 0.95,
                                     iprStatistics =
                                     matrix(twoPlNonUniNcdif[nullFocal, ], nrow =
                                            length(nullFocal)))

# Calculate cutoff point current algorithm
set.seed(29834328)

twoPlNonUniAseCurrent                <- twoPlNonUniAse
twoPlNonUniAseCurrent[['focal']]     <- twoPlNonUniAseCurrent[['focal']][nullFocal]
twoPlNonUniAseCurrent[['reference']] <- lapply(twoPlNonUniAseCurrent[['reference']][nullFocal],
                                               '*', kRatio)

twoPlNonUniIprCurrent <- Ipr(itemParameters = itemParametersNull,
                             itemCovariances = twoPlNonUniAseCurrent,
                             nReplicates = nReplicates)

twoPlNonUniNcdifCurrent <- IprNcdif(itemParameterList = twoPlNonUniIprCurrent,
                                    irtModel = "2pl", logistic = FALSE,
                                    subdivisions = 1000)
cutoffPointNonUni <- CutoffIpr(quantiles = 0.95,
                               iprStatistics = matrix(twoPlNonUniNcdifCurrent,
                                                      nrow = length(nullFocal)))

# Calculate power
powerNonUni <- data.frame(Type = "Cutoff with only focal group",
                       a = itemParameters[['focal']][, 1],
                       NCDIF = (twoPlNonUniNcdifTrue *
                                sign(itemParameters[['focal']][, 1] - 
                                     itemParametersNull[['focal']][, 1])),
                       Power = rowMeans(twoPlNonUniNcdif >
                                        cutoffPointNonUni$quantiles))

powerNonUniEach <- data.frame(Type = "Cutoff with both groups",
                           a = itemParameters[['focal']][, 1],
                           NCDIF = (twoPlNonUniNcdifTrue *
                                    sign(itemParameters[['focal']][, 1] -
                                     itemParametersNull[['focal']][, 1])),
                           Power = rowMeans(twoPlNonUniNcdif >
                                            cutoffPointEachSZNonUni$quantiles))

power80N <- which(abs(powerNonUniEach[, 'Power'] - 0.8) ==
                  min(abs(powerNonUniEach[, 'Power'] - 0.8)))

powerNonUniJ <- rbind(powerNonUni, powerNonUniEach)

powerPlotNonUniDif <- ggplot(powerNonUniJ, aes(x = a, y = Power))
powerPlotNonUniDif <- powerPlotNonUniDif + xlim(0.5, 2.0)
powerPlotNonUniDif <- powerPlotNonUniDif + geom_hline(yintercept = 0.05,
                                                linetype = 'longdash',
                                                colour = 'grey')
powerPlotNonUniDif <- powerPlotNonUniDif + geom_hline(yintercept =
                                                powerNonUniEach[power80N, 'Power'],
                                                linetype = 'longdash',
                                                colour = 'grey')
powerPlotNonUniDif <- powerPlotNonUniDif + geom_hline(yintercept =
                                                powerNonUni[power80N, 'Power'],
                                                linetype = 'longdash',
                                                colour = 'grey')
powerPlotNonUniDif <- powerPlotNonUniDif + geom_vline(xintercept =
                                                powerNonUni[power80N, 'a'],
                                                linetype = 'longdash',
                                                colour = 'grey')
powerPlotNonUniDif <- powerPlotNonUniDif + geom_vline(xintercept =
                                                -powerNonUni[power80N, 'a'],
                                                linetype = 'longdash',
                                                colour = 'grey')
powerPlotNonUniDif <- powerPlotNonUniDif + geom_line(aes(colour = Type,
                                                   linetype = Type), size = 1.0)

powerPlotNonUniNcdif <- ggplot(powerNonUniJ, aes(x = NCDIF, y = Power))
powerPlotNonUniNcdif <- powerPlotNonUniNcdif + xlim(-0.015, 0.015)
powerPlotNonUniNcdif <- powerPlotNonUniNcdif + geom_hline(yintercept = 0.05,
                                                    linetype = 'longdash',
                                                    colour = 'grey')
powerPlotNonUniNcdif <- powerPlotNonUniNcdif + geom_hline(yintercept =
                                                    powerNonUniEach[power80N, 'Power'],
                                                    linetype = 'longdash',
                                                    colour = 'grey')
powerPlotNonUniNcdif <- powerPlotNonUniNcdif + geom_hline(yintercept =
                                                    powerNonUni[power80N, 'Power'],
                                                    linetype = 'longdash',
                                                    colour = 'grey')
powerPlotNonUniNcdif <- powerPlotNonUniNcdif + geom_vline(xintercept =
                                                    powerNonUniEach[power80N, 'NCDIF'],
                                                    linetype = 'longdash',
                                                    colour = 'grey')
powerPlotNonUniNcdif <- powerPlotNonUniNcdif + geom_vline(xintercept =
                                                    -powerNonUniEach[power80N, 'NCDIF'],
                                                    linetype = 'longdash',
                                                    colour = 'grey')
powerPlotNonUniNcdif <- powerPlotNonUniNcdif + geom_line(aes(colour = Type,
                                                       linetype = Type),
                                                   size = 1.0)

powerPlotNonUniDif

powerPlotNonUniNcdif



# # Calculate bias
biasUniform <- data.frame(NCDIF = twoPlUniNcdifTrue *
                                  sign(itemParameters[['focal']][, 2]),
                          Bias = rowMeans(twoPlUniNcdif) - twoPlUniNcdifTrue)

biasUniform[, 'RelBias'] <- abs(biasUniform[, 'Bias'] / biasUniform[, 'NCDIF'])
biasUniform[, 'Type']    <- 'Uniform DIF'

biasPlotUniform <- ggplot(biasUniform, aes(x = NCDIF, y = Bias))
biasPlotUniform <- biasPlotUniform + geom_line() + ylab('Bias') + xlim(-0.015, 0.015) + ylim(0, 0.015)

biasPlotUniform

relbiasPlotUniform <- ggplot(biasUniform, aes(x = NCDIF, y = RelBias))
relbiasPlotUniform <- relbiasPlotUniform + geom_line() + ylab('Relative bias') + xlim(-0.015, 0.015)

relbiasPlotUniform
