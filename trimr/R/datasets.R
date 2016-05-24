#------------------------------------------------------------------------------
#' Example response time data set
#'
#' An example data set containing multiple participants' data for a response
#' time study involving two experimental conditions. The data set also includes
#' This is a synthetic data set and has no theoretical basis.
#'
#' @format A data frame with 20518 rows and 4 variables:
#' \describe{
#'     \item{participant}{participant identification number}
#'     \item{condition}{the experimental condition (2 in this example)}
#'     \item{rt}{response time, coded in milliseconds}
#'     \item{accuracy}{accuracy of the response; 1 = correct, 0 = error}
#'
#'     }
"exampleData"
#------------------------------------------------------------------------------


#------------------------------------------------------------------------------
#' SDs used for the recursive / moving criterion trimming methods
#'
#' A data frame containing the SDs used for each sample size as trimming
#' criterion for the nonRecursive function and the modifiedRecursive function
#'
#' @format A data frame with 97 rows and 3 columns:
#' \describe{
#'     \item{sampleSize}{Sample size of the data set being passed}
#'     \item{nonRecursive}{The standard deviation to use as the criterion for
#'     the nonRecursive function}
#'     \item{modifiedRecursive}{The standard deviation to use as the criterion
#'     for the modifiedRecursive function}
#'
#'     }
"linearInterpolation"
#------------------------------------------------------------------------------
