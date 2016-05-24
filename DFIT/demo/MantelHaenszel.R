################################################################################
# # demo.MantelHaenszel.R
# # R Versions: 2.14.1
# #
# # Author(s): Victor H. Cervantes
# #
# # Mantel Haenszel statistic assuming IRT response models
# # Description: Shows the working of the funcions in MantelHaenszel.R
# #
# # Inputs: MantelHaenszel.R
# #
# # Outputs: Example runs from each function in MantelHaenszel.R
# #
# # File history:
# #   20120412: Creation
# #   20140428: Adjustment to demo for the DFIT package
################################################################################

################################################################################
# # Load MantelHaenszel.R functions
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

exampleAbilities <- rnorm(100)

# # Dichotomous item parameters taken from Wright 2011
data(dichotomousItemParameters)
#load("../data/dichotomousItemParameters.rda")

# # Rasch model
raschParameters <- lapply(dichotomousItemParameters, function (x) matrix(x[, 2], ncol = 1))
raschParameters <- as.list(unique(as.data.frame.list(raschParameters)))
raschParameters <- lapply(raschParameters, function (x) matrix(x, ncol = 1))

# # Two parameter logistic model logistic metric
#
twoPlParameters <- lapply(dichotomousItemParameters, function (x) matrix(x[, 1:2], ncol = 2))
twoPlParameters <- unique(as.data.frame.list(twoPlParameters))
twoPlParameters <- list(focal     = as.matrix(twoPlParameters[, grep("focal", names(twoPlParameters))]),
                        reference = as.matrix(twoPlParameters[, grep("reference", names(twoPlParameters))]))
twoPlParameters <- lapply(twoPlParameters, function (x) matrix(x, ncol = 2))

# # Three parameter logistic model logistic metric
#
threePlParameters <- dichotomousItemParameters

################################################################################
# # Examples
################################################################################

mhRasch <- IrtMh(itemParameters = raschParameters, irtModel = "1pl", focalDistribution = "norm",
                 referenceDistribution = "norm", focalDistrExtra = list(mean = 0, sd = 1),
                 referenceDistrExtra = list(mean = 0, sd = 1), groupRatio = 1, logistic = TRUE)

mh2pl <- IrtMh(itemParameters = twoPlParameters, irtModel = "2pl", focalDistribution = "norm",
                 referenceDistribution = "norm", focalDistrExtra = list(mean = 0, sd = 1),
                 referenceDistrExtra = list(mean = 0, sd = 1), groupRatio = 1, logistic = FALSE)

mh3pl <- IrtMh(itemParameters = threePlParameters, irtModel = "3pl", focalDistribution = "norm",
                 referenceDistribution = "norm", focalDistrExtra = list(mean = 0, sd = 1),
                 referenceDistrExtra = list(mean = 0, sd = 1), groupRatio = 1, logistic = FALSE)

mh3plGroupImpact <- IrtMh(itemParameters = threePlParameters, irtModel = "3pl", focalDistribution = "norm",
                         referenceDistribution = "norm", focalDistrExtra = list(mean = -1, sd = 1),
                         referenceDistrExtra = list(mean = 0, sd = 1), groupRatio = 3, logistic = FALSE)


deltaRasch <- DeltaMhIrt(mhRasch)
delta2pl <- DeltaMhIrt(mh2pl)
delta3pl <- DeltaMhIrt(mh3pl)
delta3plGroupImpact <- DeltaMhIrt(mh3plGroupImpact)
