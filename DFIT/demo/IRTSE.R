###############################################################################
# # demo.IRTSE.R
# # R Versions: 2.15.0
# #
# # Author(s): Victor H. Cervantes
# #
# # Asymptotic estandard errors for item parameter estimates
# # Description: Shows the working of the funcions in IRTSE.R
# #
# # Inputs: NCDIF.R
# #
# # Outputs: Example runs from each function in IRTSE.R
# #
# # File history:
# #   20120506: Creation
# #   20140428: Adjustments to be a demo for the DFIT package
################################################################################

################################################################################
# # Load related functions
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
# # Dichotomous item parameters taken from Wright 2011
data(dichotomousItemParameters)
#load("../data/dichotomousItemParameters.rda")

# # Polytomous item parameters taken from Raju 2009
#data(polytomousItemParameters)

################################################################################
# # Dichotomous items examples
################################################################################

# # Rasch model and 1pl model
raschParameters <- dichotomousItemParameters[['focal']][, 2]
raschParameters <- unique(raschParameters)
raschParameters <- matrix(raschParameters, ncol = 1)

raschAse100 <- AseIrt(itemParameters = raschParameters, logistic = TRUE,
                      sampleSize = 100, irtModel = "1pl")

raschAse500 <- AseIrt(itemParameters = raschParameters, logistic = TRUE,
                      sampleSize = 500, irtModel = "1pl")

# # Two parameter logistic model logistic metric
#
twoPlParameters <- dichotomousItemParameters[['focal']][, -3]
twoPlParameters <- unique(twoPlParameters)
twoPlParameters <- matrix(twoPlParameters, ncol = 2)

twoPlACovNormalMetric <- AseIrt(itemParameters = twoPlParameters, logistic = FALSE,
                                sampleSize = 100, irtModel = "2pl")

twoPlACovLogisticMetric <- AseIrt(itemParameters = twoPlParameters, logistic = TRUE,
                                sampleSize = 100, irtModel = "2pl")

sqrt(twoPlACovLogisticMetric[[1]] / twoPlACovNormalMetric[[2]])

# # Three parameter logistic model logistic metric
#
threePlParameters <- dichotomousItemParameters[['focal']]

threePlACovNormalMetric <- AseIrt(itemParameters = threePlParameters, logistic = FALSE,
                                sampleSize = 100, irtModel = "3pl")

threePlACovLogisticMetric <- AseIrt(itemParameters = threePlParameters, logistic = TRUE,
                                sampleSize = 100, irtModel = "3pl")


threePlNormalMetricLi <- AseIrt(itemParameters = matrix(c(1.2, 0.5, 0.2), nrow = 1),
                                logistic = FALSE, sampleSize = 1000, irtModel = "3pl")
