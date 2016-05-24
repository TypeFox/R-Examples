#' Diabetes data.
#'
#' The diabetes data is a data frame with 11 variables and 442 measurements.
#' These are the data used in the Efron et al. (2004) paper. The data has been
#' standardized to have unit L2 norm in each column and zero mean.

#'
#'@name diabetes
#'@docType data
#'@usage diabetes
#'@format \code{diabetes} is a data frame with 11 variables (columns).
#' The first column of the data frame contains the response variable
#' (\code{diabetes$y}) which is a quantitative measure of disease progression
#' one year after baseline. The rest of the columns contain the measurements of
#' the ten explanatory variables (age, sex, body mass index, average blood
#' pressure and six blood serum registers) obtained from each of the 442
#' diabetes patients.
#'@source The orginal data are available in the lars package,
#' see http://cran.r-project.org/web/packages/lars/.
#'@references Efron, B., Hastie, T., Johnstone, I. and Tibshirani, R. (2004).
#' Least angle regression (with discussion). Annals of Statistics, 32:407--499.
#'@examples
#' library(FWDselect)
#' data(diabetes)
#' head(diabetes)
#'


NULL
