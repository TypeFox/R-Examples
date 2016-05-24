#  Import functions from non-default packages ##################################
#  Generic S3 methods must be imported in the respective file where they are 
#  used/modified so that they are included in the NAMESPACE. The following are
#  needed to pass CRAN check.
#' @importFrom graphics axis legend lines mtext par title hist plot abline
#'             arrows boxplot points text 
#' @importFrom georob georob sample.variogram
#' @importFrom stats dnorm extractAIC formula lm median model.frame cor 
#'             model.matrix residuals sd quantile var
#' @importFrom randomForest randomForest
#' @importFrom geoR practicalRange
