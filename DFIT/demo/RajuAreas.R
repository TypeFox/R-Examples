################################################################################
# # demo.RajuAreas.R
# # R Versions: 2.14.1
# #
# # Author(s): Victor H. Cervantes
# #
# # Raju's DIF Area Measures
# # Description: Shows the working of the funcions in RajuAreas.R
# #
# # Inputs: NCDIF.R
# #
# # Outputs: Example runs from each function in RajuAreas.R
# #
# # File history:
# #   20120412: Creation
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


################################################################################
# # Other global variables
################################################################################
set.seed(298374328)

# # Dichotomous item parameters taken from Wright 2011
data(dichotomousItemParameters)
#load("../data/dichotomousItemParameters.rda")

# # Polytomous item parameters taken from Raju 2009
data(polytomousItemParameters)
#load("../data/polytomousItemParameters.rda")

exampleAbilitiesSmall  <- rnorm(100)
exampleAbilitiesMedium <- rnorm(1000)
exampleAbilitiesLarge  <- rnorm(10000)


################################################################################
# # Dichotomous items examples
################################################################################

# # Rasch model and 1pl model
raschParameters <- lapply(dichotomousItemParameters, function (x) matrix(x[, 2], ncol = 1))
raschParameters <- as.list(unique(as.data.frame.list(raschParameters)))
raschParameters <- lapply(raschParameters, function (x) matrix(x, ncol = 1))

sam1pl <- SignedArea(itemParameters = raschParameters, irtModel = "1pl", 
                                subdivisions = 5000, logistic = TRUE)
uam1pl <- UnsignedArea(itemParameters = raschParameters, irtModel = "1pl",
                              subdivisions = 5000, logistic = FALSE)

# # Two parameter logistic model logistic metric
#
twoPlParameters <- lapply(dichotomousItemParameters, function (x) matrix(x[, 1:2], ncol = 2))
twoPlParameters <- unique(as.data.frame.list(twoPlParameters))
twoPlParameters <- list(focal     = as.matrix(twoPlParameters[, grep("focal", names(twoPlParameters))]),
                        reference = as.matrix(twoPlParameters[, grep("reference", names(twoPlParameters))]))
twoPlParameters <- lapply(twoPlParameters, function (x) matrix(x, ncol = 2))

# # NCDIF
sam2pl <- SignedArea(itemParameters = twoPlParameters, irtModel = "2pl",
                                subdivisions = 5000, logistic = TRUE)
uam2pl <- UnsignedArea(itemParameters = twoPlParameters, irtModel = "2pl",
                                subdivisions = 5000, logistic = TRUE)

# # Three parameter logistic model logistic metric
#
threePlParameters <- dichotomousItemParameters

sam3pl <- SignedArea(itemParameters = threePlParameters, irtModel = "3pl",
                                subdivisions = 5000, logistic = TRUE)
uam3pl <- UnsignedArea(itemParameters = threePlParameters, irtModel = "3pl",
                                subdivisions = 5000, logistic = TRUE)

################################################################################
# # Polytomous items examples
################################################################################

# # GRM

# # GRM model logistic metric
samGrm <- SignedArea(itemParameters = polytomousItemParameters, irtModel = "grm",
                                subdivisions = 5000, logistic = TRUE)
uamGrm <- UnsignedArea(itemParameters = polytomousItemParameters, irtModel = "grm",
                                subdivisions = 5000, logistic = TRUE)

# # GPCM

samGpcm <- SignedArea(itemParameters = polytomousItemParameters, irtModel = "pcm",
                                subdivisions = 5000, logistic = TRUE)
uamGpcm <- UnsignedArea(itemParameters = polytomousItemParameters, irtModel = "pcm",
                                subdivisions = 5000, logistic = TRUE)
