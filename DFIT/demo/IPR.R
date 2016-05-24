################################################################################
# # demo.IPR.R
# # R Versions: 2.15.0
# #
# # Author(s): Victor H. Cervantes
# #
# # Item parameter replication approach
# # Description: Shows the working of the funcions in IPR.R
# #
# # Outputs: Example runs from each function in IPR.R
# #
# # File history:
# #   20120403: Creation
# #   20140428: Adjustments as demo for the DFIT package
# #   20140518: Adjustments to check cutoff points
################################################################################

################################################################################
# # Load functions
################################################################################
library('DFIT')
#source("../R/IPR.R")
#source("../R/IRTSE.R")
#source("../R/MantelHaenszel.R")
#source("../R/RajuAreas.R")
#source("../R/NCDIF.R")


################################################################################
# # Other global variables
################################################################################
set.seed(29837428)

# # Dichotomous item parameters taken from Wright 2011
data(dichotomousItemParameters)
#load("../data/dichotomousItemParameters.rda")

# # Polytomous item parameters taken from Raju 2009
#data(polytomousItemParameters)

################################################################################
# # Dichotomous items examples
################################################################################

# # Rasch model
raschParameters <- lapply(dichotomousItemParameters, function (x) matrix(x[, 2], ncol = 1))
raschParameters <- as.list(unique(as.data.frame.list(raschParameters)))
raschParameters <- lapply(raschParameters, function (x) matrix(x, ncol = 1))

raschAse <- list()
raschAse[["focal"]] <- AseIrt(itemParameters = raschParameters[["focal"]], logistic = TRUE,
                      sampleSize = 500, irtModel = "1pl")
raschAse[["reference"]] <- AseIrt(itemParameters = raschParameters[["reference"]], logistic = TRUE,
                      sampleSize = 1500, irtModel = "1pl")

raschIpr      <- Ipr(itemParameters = raschParameters, itemCovariances = raschAse, nReplicates = 100)
raschNcdifIpr <- IprNcdif(itemParameterList = raschIpr, irtModel = "1pl", logistic = TRUE)
raschRUamIpr  <- IprUam(itemParameterList = raschIpr, irtModel = "1pl", logistic = TRUE)
raschMhIpr    <- IprMh(itemParameterList = raschIpr, irtModel = "1pl", logistic = TRUE)
cutoff1Pl     <- CutoffIpr(quantiles = 0.95, iprStatistics = raschNcdifIpr)

# # Two parameter logistic model
#
isEven          <- ((seq(nrow(dichotomousItemParameters[['focal']])) %% 2) == 1)
twoPlParameters <- lapply(dichotomousItemParameters, function (x) matrix(x[!isEven, 1:2], ncol = 2))
twoPlParameters <- unique(as.data.frame.list(twoPlParameters))
twoPlParameters <- list(focal     = as.matrix(twoPlParameters[, grep("focal", names(twoPlParameters))]),
                        reference = as.matrix(twoPlParameters[, grep("reference", names(twoPlParameters))]))
twoPlParameters <- lapply(twoPlParameters, function (x) matrix(x, ncol = 2))

twoPlAse <- list()
twoPlAse[["focal"]]     <- AseIrt(itemParameters = twoPlParameters[["focal"]], logistic = TRUE,
                                  sampleSize = 1000, irtModel = "2pl")
twoPlAse[["reference"]] <- AseIrt(itemParameters = twoPlParameters[["reference"]], logistic = TRUE,
                                  sampleSize = 1500, irtModel = "2pl")

twoPlIpr      <- Ipr(itemParameters = twoPlParameters, itemCovariances = twoPlAse, nReplicates = 100)
twoPlNcdifIpr <- IprNcdif(itemParameterList = twoPlIpr, irtModel = "2pl", logistic = TRUE)
twoPlRUamIpr  <- IprUam(itemParameterList = twoPlIpr, irtModel = "2pl", logistic = TRUE)
twoPlMhIpr    <- IprMh(itemParameterList = twoPlIpr, irtModel = "2pl", logistic = TRUE)
cutoff2Pl     <- CutoffIpr(quantiles = 0.95, itemParameterList = twoPlIpr, irtModel = "2pl", statistic = "uam")

# # Three parameter logistic model 
#
threePlParameters <- dichotomousItemParameters
isNot3Pl          <- ((dichotomousItemParameters[['focal']][, 3] == 0) | (dichotomousItemParameters[['reference']][, 3] == 0))

threePlParameters[['focal']]          <- threePlParameters[['focal']][!isNot3Pl, ]
threePlParameters[['reference']]      <- threePlParameters[['reference']][!isNot3Pl, ]
threePlParameters[['focal']][, 3]     <- threePlParameters[['focal']][, 3] + 0.1
threePlParameters[['reference']][, 3] <- threePlParameters[['reference']][, 3] + 0.1
threePlParameters[['focal']][, 2]     <- threePlParameters[['focal']][, 2] + 1.5
threePlParameters[['reference']][, 2] <- threePlParameters[['reference']][, 2] + 1.5
threePlParameters[['focal']]          <- threePlParameters[['focal']][-c(12, 16, 28), ]
threePlParameters[['reference']]      <- threePlParameters[['reference']][-c(12, 16, 28), ]


threePlAse <- list()
threePlAse[["focal"]]     <- AseIrt(itemParameters = threePlParameters[["focal"]], logistic = TRUE,
                                    sampleSize = 10000, irtModel = "3pl")
threePlAse[["reference"]] <- AseIrt(itemParameters = threePlParameters[["reference"]], logistic = TRUE,
                                    sampleSize = 15000, irtModel = "3pl")

set.seed(41568)
threePlIpr <- Ipr(itemParameters = threePlParameters, itemCovariances = threePlAse, nReplicates = 100)

threePlIpr <- Bound3PlIpr(threePlIpr)

threePlNcdifIpr <- IprNcdif(itemParameterList = threePlIpr, irtModel = "3pl", logistic = TRUE)
threePlRUamIpr  <- IprUam(itemParameterList = threePlIpr, irtModel = "3pl", logistic = TRUE)
threePlMhIpr    <- IprMh(itemParameterList = threePlIpr, irtModel = "3pl", logistic = TRUE)

set.seed(41568)
cutoff3Pl       <- CutoffIpr(quantiles = 0.95, itemParameters = threePlParameters, itemCovariances = threePlAse,
                             irtModel = "3pl", statistic = "mh", nReplicates = 100)

