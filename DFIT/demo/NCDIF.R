################################################################################
# # demo.NCDIF.R
# # R Versions: 2.14.1
# #
# # Author(s): Victor H. Cervantes
# #
# # General NCDIF function
# # Description: Shows the working of the funcions in NCDIF.R
# #
# # Inputs: NCDIF.R
# #
# # Outputs: Example runs from each function in NCDIF.R
# #
# # File history:
# #   20120403: Creation
# #   20140428: Adjustments to be a demo for the DFIT package
################################################################################

################################################################################
# # Load NCDIF.R functions
################################################################################
library("DFIT")
#source("../R/IPR.R")
#source("../R/IRTSE.R")
#source("../R/MantelHaenszel.R")
#source("../R/RajuAreas.R")
#source("../R/NCDIF.R")
#library(ggplot2)


################################################################################
# # Other global variables
################################################################################
set.seed(298374328)

exampleAbilitiesSmall  <- rnorm(100)
exampleAbilitiesMedium <- rnorm(1000)
exampleAbilitiesLarge  <- rnorm(10000)

# # Dichotomous item parameters taken from Wright 2011
data(dichotomousItemParameters)
#load("../data/dichotomousItemParameters.rda")

# # Polytomous item parameters taken from Raju 2009
data(polytomousItemParameters)
#load("../data/polytomousItemParameters.rda")

################################################################################
# # Dichotomous items examples
################################################################################

# # Rasch model
raschParameters <- lapply(dichotomousItemParameters, function (x) matrix(x[, 2], ncol = 1))
raschParameters <- as.list(unique(as.data.frame.list(raschParameters)))
raschParameters <- lapply(raschParameters, function (x) matrix(x, ncol = 1))

# # NCDIF
ncdifRaschTheoric      <- Ncdif(itemParameters = raschParameters, irtModel = "1pl", focalAbilities = NULL,
                                focalDistribution = "norm", subdivisions = 5000, logistic = TRUE)
ncdifRaschSmallSample  <- Ncdif(itemParameters = raschParameters, irtModel = "1pl", focalAbilities = exampleAbilitiesSmall,
                                focalDistribution = "norm", subdivisions = 5000, logistic = TRUE)
ncdifRaschMediumSample <- Ncdif(itemParameters = raschParameters, irtModel = "1pl",
                                focalAbilities = exampleAbilitiesMedium,
                                focalDistribution = "norm", subdivisions = 5000, logistic = TRUE)
ncdifRaschLargeSample  <- Ncdif(itemParameters = raschParameters, irtModel = "1pl",
                                focalAbilities = exampleAbilitiesLarge,
                                focalDistribution = "norm", subdivisions = 5000, logistic = TRUE)

# # CDIF and DFT
cdifRaschTheoric           <- Cdif(itemParameters = raschParameters, irtModel = "1pl", focalAbilities = NULL,
                                   focalDistribution = "norm", subdivisions = 5000, logistic = TRUE)
dtfRaschTheoricWithCdif    <- Dtf(cdif = cdifRaschTheoric)
dtfRaschTheoricWithoutCdif <- Dtf(cdif = NULL, itemParameters = raschParameters, irtModel = "1pl", focalAbilities = NULL,
                                  focalDistribution = "norm", subdivisions = 5000, logistic = TRUE)

# # One parameter logistic model normal metric
# # NCDIF
ncdif1plTheoric      <- Ncdif(itemParameters = raschParameters, irtModel = "1pl", focalAbilities = NULL,
                              focalDistribution = "norm", subdivisions = 5000, logistic = FALSE)
ncdif1plSmallSample  <- Ncdif(itemParameters = raschParameters, irtModel = "1pl", focalAbilities = exampleAbilitiesSmall,
                              focalDistribution = "norm", subdivisions = 5000, logistic = FALSE)
ncdif1plMediumSample <- Ncdif(itemParameters = raschParameters, irtModel = "1pl",
                              focalAbilities = exampleAbilitiesMedium,
                              focalDistribution = "norm", subdivisions = 5000, logistic = FALSE)
ncdif1plLargeSample  <- Ncdif(itemParameters = raschParameters, irtModel = "1pl",
                              focalAbilities = exampleAbilitiesLarge,
                              focalDistribution = "norm", subdivisions = 5000, logistic = FALSE)

# # CDIF and DFT
cdif1plNormMetricTheoric        <- Cdif(itemParameters = raschParameters, irtModel = "1pl", focalAbilities = NULL,
                                        focalDistribution = "norm", subdivisions = 5000, logistic = FALSE)
cdif1plNormMetricLargeSample    <- Cdif(itemParameters = raschParameters, irtModel = "1pl",
                                        focalAbilities = exampleAbilitiesLarge,
                                        focalDistribution = "norm", subdivisions = 5000, logistic = FALSE)
dtf1plNormMetricTheoricWithCdif <- Dtf(cdif = cdif1plNormMetricTheoric)
dtf1plNormMetricWithoutCdif     <- Dtf(cdif = NULL, itemParameters = raschParameters, irtModel = "1pl", focalAbilities = NULL,
                                       focalDistribution = "norm", subdivisions = 5000, logistic = FALSE)
dtf1plNormMetricLargeSample     <- Dtf(cdif = cdif1plNormMetricLargeSample)


# # Plot examples for item with and without uniform DIF under 1pl model

# # No DIF item
PlotNcdif(iiItem = 5, itemParameters = raschParameters, irtModel = "1pl", plotDensity = TRUE,
          focalDensityText = "Focal group density (theoretical)",
          main = "Item 5. NO DIF. Rasch model")

# # DIF item close to focal distribution
PlotNcdif(iiItem = 7, itemParameters = raschParameters, irtModel = "1pl", plotDensity = TRUE,
          main = "Item 7. Uniform DIF. Rasch model")
PlotNcdif(iiItem = 7, itemParameters = raschParameters, irtModel = "1pl", plotDensity = FALSE,
          main = "Item 7. Uniform DIF. Rasch model")

# # DIF item far from focal distribution
PlotNcdif(iiItem = 3, itemParameters = raschParameters, irtModel = "1pl", plotDensity = TRUE)
PlotNcdif(iiItem = 3, itemParameters = raschParameters, irtModel = "1pl", plotDensity = FALSE)



# # Two parameter logistic model logistic metric
#
twoPlParameters <- lapply(dichotomousItemParameters, function (x) matrix(x[, 1:2], ncol = 2))
twoPlParameters <- unique(as.data.frame.list(twoPlParameters))
twoPlParameters <- list(focal     = as.matrix(twoPlParameters[, grep("focal", names(twoPlParameters))]),
                        reference = as.matrix(twoPlParameters[, grep("reference", names(twoPlParameters))]))
twoPlParameters <- lapply(twoPlParameters, function (x) matrix(x, ncol = 2))

# # NCDIF
ncdif2plLogMetricTheoric      <- Ncdif(itemParameters = twoPlParameters, irtModel = "2pl", focalAbilities = NULL,
                                       focalDistribution = "norm", subdivisions = 5000, logistic = TRUE)
ncdif2plLogMetricSmallSample  <- Ncdif(itemParameters = twoPlParameters, irtModel = "2pl", focalAbilities = exampleAbilitiesSmall,
                                       focalDistribution = "norm", subdivisions = 5000, logistic = TRUE)
ncdif2plLogMetricMediumSample <- Ncdif(itemParameters = twoPlParameters, irtModel = "2pl",
                                       focalAbilities = exampleAbilitiesMedium,
                                       focalDistribution = "norm", subdivisions = 5000, logistic = TRUE)
ncdif2plLogMetricLargeSample  <- Ncdif(itemParameters = twoPlParameters, irtModel = "2pl",
                                       focalAbilities = exampleAbilitiesLarge,
                                       focalDistribution = "norm", subdivisions = 5000, logistic = TRUE)

# # CDIF and DFT
cdif2plLogMetricTheoric        <- Cdif(itemParameters = twoPlParameters, irtModel = "2pl", focalAbilities = NULL,
                                       focalDistribution = "norm", subdivisions = 5000, logistic = TRUE)
cdif2plLogMetricLargeSample    <- Cdif(itemParameters = twoPlParameters, irtModel = "2pl",
                                       focalAbilities = exampleAbilitiesLarge,
                                       focalDistribution = "norm", subdivisions = 5000, logistic = TRUE)
dtf2plLogMetricTheoricWithCdif <- Dtf(cdif = cdif2plLogMetricTheoric)
dtf2plLogMetricWithoutCdif     <- Dtf(cdif = NULL, itemParameters = twoPlParameters, irtModel = "2pl", focalAbilities = NULL,
                                      focalDistribution = "norm", subdivisions = 5000, logistic = TRUE)
dtf2plLogMetricLargeSample     <- Dtf(cdif = cdif2plLogMetricLargeSample)

# # Two parameter logistic model normal metric
# # NCDIF
ncdif2plNormMetricTheoric      <- Ncdif(itemParameters = twoPlParameters, irtModel = "2pl", focalAbilities = NULL,
                                        focalDistribution = "norm", subdivisions = 5000, logistic = FALSE)
ncdif2plNormMetricSmallSample  <- Ncdif(itemParameters = twoPlParameters, irtModel = "2pl", focalAbilities = exampleAbilitiesSmall,
                                        focalDistribution = "norm", subdivisions = 5000, logistic = FALSE)
ncdif2plNormMetricMediumSample <- Ncdif(itemParameters = twoPlParameters, irtModel = "2pl",
                                        focalAbilities = exampleAbilitiesMedium,
                                        focalDistribution = "norm", subdivisions = 5000, logistic = FALSE)
ncdif2plNormMetricLargeSample  <- Ncdif(itemParameters = twoPlParameters, irtModel = "2pl",
                                        focalAbilities = exampleAbilitiesLarge,
                                        focalDistribution = "norm", subdivisions = 5000, logistic = FALSE)

# # CDIF and DFT
cdif2plNormMetricTheoric        <- Cdif(itemParameters = twoPlParameters, irtModel = "2pl", focalAbilities = NULL,
                                        focalDistribution = "norm", subdivisions = 5000, logistic = FALSE)
cdif2plNormMetricLargeSample    <- Cdif(itemParameters = twoPlParameters, irtModel = "2pl",
                                        focalAbilities = exampleAbilitiesLarge,
                                        focalDistribution = "norm", subdivisions = 5000, logistic = FALSE)
dtf2plNormMetricTheoricWithCdif <- Dtf(cdif = cdif2plNormMetricTheoric)
dtf2plNormMetricWithoutCdif     <- Dtf(cdif = NULL, itemParameters = twoPlParameters, irtModel = "2pl", focalAbilities = NULL,
                                       focalDistribution = "norm", subdivisions = 5000, logistic = FALSE)
dtf2plNormMetricLargeSample     <- Dtf(cdif = cdif2plNormMetricLargeSample)

# # Plot examples for item with non uniform and mixed DIF under 2pl model

# # Non Uniform DIF item close to focal distribution
PlotNcdif(iiItem = 19, itemParameters = twoPlParameters, irtModel = "2pl", plotDensity = TRUE,
          main = "Item 19. Non uniform DIF. 2PL")
PlotNcdif(iiItem = 19, itemParameters = twoPlParameters, irtModel = "2pl", plotDensity = FALSE,
          main = "Item 19. Non uniform DIF. 2PL")

# # Non Uniform DIF item far from focal distribution
PlotNcdif(iiItem = 2, itemParameters = twoPlParameters, irtModel = "2pl", plotDensity = FALSE,
          main = "Item 2. Non uniform DIF. 2PL")

# # Mixed DIF item close to focal distribution
PlotNcdif(iiItem = 15, itemParameters = twoPlParameters, irtModel = "2pl", plotDensity = FALSE,
          main = "Item 19. Non uniform DIF. 2PL")

# # Mixed DIF item far from focal distribution
PlotNcdif(iiItem = 23, itemParameters = twoPlParameters, irtModel = "2pl", plotDensity = FALSE,
          main = "Item 2. Non uniform DIF. 2PL")



# # Three parameter logistic model logistic metric
#
threePlParameters <- dichotomousItemParameters

# # NCDIF

ncdif3plLogMetricTheoric      <- Ncdif(itemParameters = threePlParameters, irtModel = "3pl", focalAbilities = NULL,
                                       focalDistribution = "norm", subdivisions = 5000, logistic = TRUE)
ncdif3plLogMetricSmallSample  <- Ncdif(itemParameters = threePlParameters, irtModel = "3pl", focalAbilities = exampleAbilitiesSmall,
                                       focalDistribution = "norm", subdivisions = 5000, logistic = TRUE)
ncdif3plLogMetricMediumSample <- Ncdif(itemParameters = threePlParameters, irtModel = "3pl",
                                       focalAbilities = exampleAbilitiesMedium,
                                       focalDistribution = "norm", subdivisions = 5000, logistic = TRUE)
ncdif3plLogMetricLargeSample  <- Ncdif(itemParameters = threePlParameters, irtModel = "3pl",
                                       focalAbilities = exampleAbilitiesLarge,
                                       focalDistribution = "norm", subdivisions = 5000, logistic = TRUE)

# # CDIF and DFT
cdif3plLogMetricTheoric        <- Cdif(itemParameters = threePlParameters, irtModel = "3pl", focalAbilities = NULL,
                                       focalDistribution = "norm", subdivisions = 5000, logistic = TRUE)
cdif3plLogMetricLargeSample    <- Cdif(itemParameters = threePlParameters, irtModel = "3pl",
                                       focalAbilities = exampleAbilitiesLarge,
                                       focalDistribution = "norm", subdivisions = 5000, logistic = TRUE)
dtf3plLogMetricTheoricWithCdif <- Dtf(cdif = cdif3plLogMetricTheoric)
dtf3plLogMetricWithoutCdif     <- Dtf(cdif = NULL, itemParameters = threePlParameters, irtModel = "3pl", focalAbilities = NULL,
                                      focalDistribution = "norm", subdivisions = 5000, logistic = TRUE)
dtf3plLogMetricLargeSample     <- Dtf(cdif = cdif3plLogMetricLargeSample)

# # Three parameter logistic model normal metric
# # NCDIF
ncdif3plNormMetricTheoric      <- Ncdif(itemParameters = threePlParameters, irtModel = "3pl", focalAbilities = NULL,
                                        focalDistribution = "norm", subdivisions = 5000, logistic = FALSE)
ncdif3plNormMetricSmallSample  <- Ncdif(itemParameters = threePlParameters, irtModel = "3pl", focalAbilities = exampleAbilitiesSmall,
                                        focalDistribution = "norm", subdivisions = 5000, logistic = FALSE)
ncdif3plNormMetricMediumSample <- Ncdif(itemParameters = threePlParameters, irtModel = "3pl",
                                        focalAbilities = exampleAbilitiesMedium,
                                        focalDistribution = "norm", subdivisions = 5000, logistic = FALSE)
ncdif3plNormMetricLargeSample  <- Ncdif(itemParameters = threePlParameters, irtModel = "3pl",
                                        focalAbilities = exampleAbilitiesLarge,
                                        focalDistribution = "norm", subdivisions = 5000, logistic = FALSE)

# # CDIF and DFT
cdif3plNormMetricTheoric        <- Cdif(itemParameters = threePlParameters, irtModel = "3pl", focalAbilities = NULL,
                                        focalDistribution = "norm", subdivisions = 5000, logistic = FALSE)
cdif3plNormMetricLargeSample    <- Cdif(itemParameters = threePlParameters, irtModel = "3pl",
                                        focalAbilities = exampleAbilitiesLarge,
                                        focalDistribution = "norm", subdivisions = 5000, logistic = FALSE)
dtf3plNormMetricTheoricWithCdif <- Dtf(cdif = cdif3plNormMetricTheoric)
dtf3plNormMetricWithoutCdif     <- Dtf(cdif = NULL, itemParameters = threePlParameters, irtModel = "3pl", focalAbilities = NULL,
                                       focalDistribution = "norm", subdivisions = 5000, logistic = FALSE)
dtf3plNormMetricLargeSample     <- Dtf(cdif = cdif3plNormMetricLargeSample)

# # Plot examples for item with non uniform and different guessing paramters DIF under 3pl model

# # Non Uniform + != guess DIF item close to focal distribution
PlotNcdif(iiItem = 30, itemParameters = threePlParameters, irtModel = "3pl", plotDensity = FALSE,
          main = "Item 30 Non uniform and different guessing DIF. 3PL")

# # != guess DIF item close to focal distribution
PlotNcdif(iiItem = 25, itemParameters = threePlParameters, irtModel = "3pl", plotDensity = FALSE,
          main = "Item 25 Different guessing DIF. 3PL")

# # Non Uniform + != guess DIF item far from focal distribution
PlotNcdif(iiItem = 22, itemParameters = threePlParameters, irtModel = "3pl", plotDensity = FALSE,
          main = "Item 22. Non uniform and different guessing DIF. 3PL")

# # Non Uniform - != guess DIF item close to focal distribution
PlotNcdif(iiItem = 46, itemParameters = threePlParameters, irtModel = "3pl", plotDensity = FALSE,
          main = "Item 46 Non uniform and different guessing DIF. 3PL")

# # Non Uniform - != guess DIF item far from focal distribution
PlotNcdif(iiItem = 38, itemParameters = threePlParameters, irtModel = "3pl", plotDensity = FALSE,
          main = "Item 38 Non uniform and different guessing DIF. 3PL")

################################################################################
# # Polytomous items examples
################################################################################

# # GRM

# # GRM model logistic metric
ncdifGrmLogMetricTheoric      <- Ncdif(itemParameters = polytomousItemParameters, irtModel = "grm", focalAbilities = NULL,
                                       focalDistribution = "norm", subdivisions = 5000, logistic = TRUE)
ncdifGrmLogMetricSmallSample  <- Ncdif(itemParameters = polytomousItemParameters, irtModel = "grm", focalAbilities = exampleAbilitiesSmall,
                                       focalDistribution = "norm", subdivisions = 5000, logistic = TRUE)
ncdifGrmLogMetricMediumSample <- Ncdif(itemParameters = polytomousItemParameters, irtModel = "grm",
                                       focalAbilities = exampleAbilitiesMedium,
                                       focalDistribution = "norm", subdivisions = 5000, logistic = TRUE)
ncdifGrmLogMetricLargeSample  <- Ncdif(itemParameters = polytomousItemParameters, irtModel = "grm",
                                       focalAbilities = exampleAbilitiesLarge,
                                       focalDistribution = "norm", subdivisions = 5000, logistic = TRUE)

# # GRM parameter logistic model normal metric
ncdifGrmNormMetricTheoric      <- Ncdif(itemParameters = polytomousItemParameters, irtModel = "grm", focalAbilities = NULL,
                                        focalDistribution = "norm", subdivisions = 5000, logistic = FALSE)
ncdifGrmNormMetricSmallSample  <- Ncdif(itemParameters = polytomousItemParameters, irtModel = "grm", focalAbilities = exampleAbilitiesSmall,
                                        focalDistribution = "norm", subdivisions = 5000, logistic = FALSE)
ncdifGrmNormMetricMediumSample <- Ncdif(itemParameters = polytomousItemParameters, irtModel = "grm",
                                        focalAbilities = exampleAbilitiesMedium,
                                        focalDistribution = "norm", subdivisions = 5000, logistic = FALSE)
ncdifGrmNormMetricLargeSample  <- Ncdif(itemParameters = polytomousItemParameters, irtModel = "grm",
                                        focalAbilities = exampleAbilitiesLarge,
                                        focalDistribution = "norm", subdivisions = 5000, logistic = FALSE)

# # CDIF and DFT
cdifGrmNormMetricTheoric        <- Cdif(itemParameters = polytomousItemParameters, irtModel = "grm", focalAbilities = NULL,
                                        focalDistribution = "norm", subdivisions = 5000, logistic = FALSE)
cdifGrmNormMetricLargeSample    <- Cdif(itemParameters = polytomousItemParameters, irtModel = "grm",
                                        focalAbilities = exampleAbilitiesLarge,
                                        focalDistribution = "norm", subdivisions = 5000, logistic = FALSE)
dtfGrmNormMetricTheoricWithCdif <- Dtf(cdif = cdifGrmNormMetricTheoric)
dtfGrmNormMetricWithoutCdif     <- Dtf(cdif = NULL, itemParameters = polytomousItemParameters, irtModel = "grm", focalAbilities = NULL,
                                       focalDistribution = "norm", subdivisions = 5000, logistic = FALSE)
dtfGrmNormMetricLargeSample     <- Dtf(cdif = cdifGrmNormMetricLargeSample)


# # Plot examples for item with uniform, non uniform and mixed DIF under GRM model

# # Uniform DIF item
PlotNcdif(iiItem = 7, itemParameters = polytomousItemParameters, irtModel = "grm", plotDensity = FALSE,
          main = "Item 7. Uniform DIF. GRM", ylab = "Expected score")

# # Non Uniform DIF item
PlotNcdif(iiItem = 10, itemParameters = polytomousItemParameters, irtModel = "grm", plotDensity = FALSE,
          main = "Item 10. Non uniform DIF. GRM", ylab = "Expected score")

# # Mixed DIF item
PlotNcdif(iiItem = 4, itemParameters = polytomousItemParameters, irtModel = "grm", plotDensity = FALSE,
          main = "Item 4. Mixed DIF. GRM", ylab = "Expected score")


# # GPCM

# # GPCM model logistic metric
ncdifGpcmLogMetricTheoric      <- Ncdif(itemParameters = polytomousItemParameters, irtModel = "pcm", focalAbilities = NULL,
                                        focalDistribution = "norm", subdivisions = 5000, logistic = TRUE)
ncdifGpcmLogMetricSmallSample  <- Ncdif(itemParameters = polytomousItemParameters, irtModel = "pcm", focalAbilities = exampleAbilitiesSmall,
                                        focalDistribution = "norm", subdivisions = 5000, logistic = TRUE)
ncdifGpcmLogMetricMediumSample <- Ncdif(itemParameters = polytomousItemParameters, irtModel = "pcm",
                                        focalAbilities = exampleAbilitiesMedium,
                                        focalDistribution = "norm", subdivisions = 5000, logistic = TRUE)
ncdifGpcmLogMetricLargeSample  <- Ncdif(itemParameters = polytomousItemParameters, irtModel = "pcm",
                                        focalAbilities = exampleAbilitiesLarge,
                                        focalDistribution = "norm", subdivisions = 5000, logistic = TRUE)

# # GPCM parameter logistic model normal metric
ncdifGpcmNormMetricTheoric      <- Ncdif(itemParameters = polytomousItemParameters, irtModel = "pcm", focalAbilities = NULL,
                                         focalDistribution = "norm", subdivisions = 5000, logistic = FALSE)
ncdifGpcmNormMetricSmallSample  <- Ncdif(itemParameters = polytomousItemParameters, irtModel = "pcm", focalAbilities = exampleAbilitiesSmall,
                                         focalDistribution = "norm", subdivisions = 5000, logistic = FALSE)
ncdifGpcmNormMetricMediumSample <- Ncdif(itemParameters = polytomousItemParameters, irtModel = "pcm",
                                         focalAbilities = exampleAbilitiesMedium,
                                         focalDistribution = "norm", subdivisions = 5000, logistic = FALSE)
ncdifGpcmNormMetricLargeSample  <- Ncdif(itemParameters = polytomousItemParameters, irtModel = "pcm",
                                         focalAbilities = exampleAbilitiesLarge,
                                         focalDistribution = "norm", subdivisions = 5000, logistic = FALSE)

# # CDIF and DFT
cdifGpcmNormMetricTheoric        <- Cdif(itemParameters = polytomousItemParameters, irtModel = "pcm", focalAbilities = NULL,
                                        focalDistribution = "norm", subdivisions = 5000, logistic = FALSE)
cdifGpcmNormMetricLargeSample    <- Cdif(itemParameters = polytomousItemParameters, irtModel = "pcm",
                                        focalAbilities = exampleAbilitiesLarge,
                                        focalDistribution = "norm", subdivisions = 5000, logistic = FALSE)
dtfGpcmNormMetricTheoricWithCdif <- Dtf(cdif = cdifGpcmNormMetricTheoric)
dtfGpcmNormMetricWithoutCdif     <- Dtf(cdif = NULL, itemParameters = polytomousItemParameters, irtModel = "pcm", focalAbilities = NULL,
                                       focalDistribution = "norm", subdivisions = 5000, logistic = FALSE)
dtfGpcmNormMetricLargeSample     <- Dtf(cdif = cdifGpcmNormMetricLargeSample)

# # Plot examples for item with uniform, non uniform and mixed DIF under GPCM model

# # Uniform DIF item
PlotNcdif(iiItem = 7, itemParameters = polytomousItemParameters, irtModel = "pcm", plotDensity = FALSE,
          main = "Item 7. Uniform DIF. GPCM", ylab = "Expected score")

# # Non Uniform DIF item
PlotNcdif(iiItem = 10, itemParameters = polytomousItemParameters, irtModel = "pcm", plotDensity = FALSE,
          main = "Item 10. Non uniform DIF. GPCM", ylab = "Expected score")

# # Mixed DIF item
PlotNcdif(iiItem = 4, itemParameters = polytomousItemParameters, irtModel = "pcm", plotDensity = FALSE,
          main = "Item 4. Mixed DIF. GPCM", ylab = "Expected score")

