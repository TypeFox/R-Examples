#' @name supermarkets
#' @title Dutch Supermarkets Data Set
#' @description This data set relates to 220 consumers rating 10 Dutch supermarket chains according 
#' to 8 variables. A rating scale from 1 to 10 was used.
#' @docType data
#' @usage supermarkets
#' @format A three-way array with supermarkets in the first dimension, variables in the second and 
#' consumers in the third dimension.
#' @source Michel van de Velden
#' @examples
#' data("supermarkets")
#' fit <- lsbclust(data = supermarkets, nclust = 6, fixed = "rows", nstart = 2)
NULL