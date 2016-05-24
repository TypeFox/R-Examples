#' @name class.probabilities.interaction.test
#' @title Class probabilities for Q-matrix containing interaction between attributes
#' @description This data set list the class probabities for Q-matrix containing interaction between attributes.
#' This data set is used to check output of \code{ScoreDCM}
#' @docType data
#' @usage class.probabilities.interaction.test
#' @format a data frame containing 15 students and 3 attributes. 1 student per row
#' @author Margi Dubal \email{margidubal@@gmail.com} & Diane Losardo \email{dlosardo@@amplify.com}
NULL

#' @name class.probabilities.test
#' @title Class probabilities for Q-matrix containing no interaction between attributes
#' @description This data set list the class probabities for Q-matrix containing no interaction between attributes.
#' This data set is used to check output of \code{ScoreDCM}
#' @docType data
#' @usage class.probabilities.test
#' @format a data frame containing 15 students and 3 attributes. 1 student per row
#' @author Margi Dubal \email{margidubal@@gmail.com} & Diane Losardo \email{dlosardo@@amplify.com}
NULL

#' @name observations.test
#' @title Obervations
#' @description This data set list the dichotomous responses of students for a test containing 11 items.
#' @docType data
#' @usage observations.test
#' @format a data frame containing 15 students and 11 items. 1 student per row
#' @author Margi Dubal \email{margidubal@@gmail.com} & Diane Losardo \email{dlosardo@@amplify.com}
NULL

#' @name qmatrix.test
#' @title Q-matrix
#' @description This data set list 1s and 0s indicating relation between items and attributes. This matrix specifies which items are required
#'  for mastery of each attribute (i.e., latent variable). The matrix must be a size of \code{nitems X nattributes}
#' @docType data
#' @usage qmatrix.test
#' @format a matrix containing 11 items and 3 attributes
#' @author Margi Dubal \email{margidubal@@gmail.com} & Diane Losardo \email{dlosardo@@amplify.com}
NULL

#' @name qmatrix.test.interaction
#' @title Q-matrix with interaction
#' @description This data set list 1s and 0s indicating relation between items and attributes. This matrix specifies which items are required
#'  for mastery of each attribute (i.e., latent variable). The matrix must be a size of \code{nitems X nattributes}
#' @docType data
#' @usage qmatrix.interaction.test
#' @format a matrix containing 11 items and 3 attributes
#' @author Margi Dubal \email{margidubal@@gmail.com} & Diane Losardo \email{dlosardo@@amplify.com}
NULL

#' @name parameter.means.DCM.Mplus.test
#' @title Parameter estimates calibrated using Mplus for fully saturated model.
#' @description a numerical vector of calibrated item and structural parameters. Values must be in the order of \code{\link{GetParameterNames}} if
#' parametrization method is Mplus and non-kernel parameters are used. If kernel parameters values are used must be in order of \code{\link{GetKernelParameterNames}}
#' @docType data
#' @usage parameter.means.DCM.Mplus.test
#' @format a vector or a dataframe
#' @author Margi Dubal \email{margidubal@@gmail.com} & Diane Losardo \email{dlosardo@@amplify.com}
NULL

#' @name parameter.means.DCM.Mplus.interaction.test
#' @title Parameter estimates calibrated using Mplus for fully saturated model and Q-matrix with interaction.
#' @description a numerical vector of calibrated item and structural parameters. Values must be in the order of \code{\link{GetParameterNames}} if
#' parametrization method is Mplus and non-kernel parameters are used. If kernel parameters values are used must be in order of \code{\link{GetKernelParameterNames}}
#' @docType data
#' @usage parameter.means.DCM.Mplus.interaction.test
#' @format a vector or a dataframe
#' @author Margi Dubal \email{margidubal@@gmail.com} & Diane Losardo \email{dlosardo@@amplify.com}
NULL

#' @name parameter.acov.DCM.Mplus.test
#' @title Covariance matrix of parameter estimates calibrated using Mplus for fully saturated model.
#' @description  matrix of covariances of all model parameters. If \code{NULL} (the default) model parameters are not randomized 
#' for each iteration of MCMC
#' @docType data
#' @usage parameter.acov.DCM.Mplus.test
#' @format A data frame
#' @author Margi Dubal \email{margidubal@@gmail.com} & Diane Losardo \email{dlosardo@@amplify.com}
NULL

#' @name parameter.acov.DCM.Mplus.interaction.test
#' @title Covariance matrix of parameter estimates calibrated using Mplus for fully saturated model and Q-matrix with interaction.
#' @description  matrix of covariances of all model parameters. If \code{NULL} (the default) model parameters are not randomized 
#' for each iteration of MCMC
#' @docType data
#' @usage parameter.acov.DCM.Mplus.interaction.test
#' @format A data frame
#' @author Margi Dubal \email{margidubal@@gmail.com} & Diane Losardo \email{dlosardo@@amplify.com}
NULL

#' @name parameter.means.names.DCM.Mplus.test
#' @title Names of parameter estimates calibrated using Mplus for fully saturated model and Q-matrix with  no interaction.
#' @description  a vector of string containing names of calibrated parameter estimates.
#' @docType data
#' @usage parameter.means.names.DCM.Mplus.test
#' @format A vector of string
#' @author Margi Dubal \email{margidubal@@gmail.com} & Diane Losardo \email{dlosardo@@amplify.com}
NULL

#' @name parameter.means.names.DCM.Mplus.interaction.test
#' @title Names of parameter estimates calibrated using Mplus for fully saturated model and Q-matrix with interaction.
#' @description  a vector of string containing names of calibrated parameter estimates.
#' @docType data
#' @usage parameter.means.names.DCM.Mplus.interaction.test
#' @format A vector of string
#' @author Margi Dubal \email{margidubal@@gmail.com} & Diane Losardo \email{dlosardo@@amplify.com}
NULL


#' @name parameter.means.DCM.kernel.Mplus.test
#' @title Kernel parameter estimates calibrated using Mplus for fully saturated model.
#' @description a numerical vector of calibrated item and structural kernel parameters. Values must be in the order of \code{\link{GetParameterNames}} if
#' parametrization method is Mplus and non-kernel parameters are used. If kernel parameters values are used must be in order of \code{\link{GetKernelParameterNames}}
#' @docType data
#' @usage parameter.means.DCM.kernel.Mplus.test
#' @format a vector or a dataframe
#' @author Margi Dubal \email{margidubal@@gmail.com} & Diane Losardo \email{dlosardo@@amplify.com}
NULL

#' @name parameter.means.DCM.kernel.Mplus.interaction.test
#' @title Kernel parameter estimates calibrated using Mplus for fully saturated model and Q-matrix with interaction.
#' @description a numerical vector of calibrated item and structural parameters. Values must be in the order of \code{\link{GetParameterNames}} if
#' parametrization method is Mplus and non-kernel parameters are used. If kernel parameters values are used must be in order of \code{\link{GetKernelParameterNames}}
#' @docType data
#' @usage parameter.means.DCM.kernel.Mplus.interaction.test
#' @format a vector or a dataframe
#' @author Margi Dubal \email{margidubal@@gmail.com} & Diane Losardo \email{dlosardo@@amplify.com}
NULL

#' @name parameter.means.names.DCM.kernel.Mplus.interaction.test
#' @title Names of kernel parameter estimates calibrated using Mplus for fully saturated model and Q-matrix with interaction.
#' @description  a vector of string containing names of calibrated parameter estimates.
#' @docType data
#' @usage parameter.means.names.DCM.kernel.Mplus.interaction.test
#' @format A vector of string
#' @author Margi Dubal \email{margidubal@@gmail.com} & Diane Losardo \email{dlosardo@@amplify.com}
NULL

#' @name parameter.means.NCRUM.interaction.test
#' @title Parameter estimates calibrated using Mplus for NC-RUM model and Q-matrix with interaction.
#' @description a numerical vector of calibrated item and structural parameters. Values must be in the order of \code{\link{GetParameterNames}} if
#' parametrization method is Mplus and non-kernel parameters are used. If kernel parameters values are used must be in order of \code{\link{GetKernelParameterNames}}
#' @docType data
#' @usage parameter.means.NCRUM.interaction.test
#' @format a vector or a dataframe
#' @author Margi Dubal \email{margidubal@@gmail.com} & Diane Losardo \email{dlosardo@@amplify.com}
NULL


#' @name parameter.means.names.NCRUM.interaction.test
#' @title Names of parameter estimates calibrated using Mplus for fully saturated model and Q-matrix with interaction.
#' @description  a vector of string containing names of calibrated parameter estimates.
#' @docType data
#' @usage parameter.means.names.NCRUM.interaction.test
#' @format A vector of string
#' @author Margi Dubal \email{margidubal@@gmail.com} & Diane Losardo \email{dlosardo@@amplify.com}
NULL