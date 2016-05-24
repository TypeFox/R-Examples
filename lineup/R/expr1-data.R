#' Example gene expression data
#'
#' Matrices of simulated gene expression data, each for 98 individuals
#' at 5,000 genes. Think of \code{expr1} and \code{expr2} as
#' expression data on two different tissues.
#'
#' @name expr-data
#' @aliases expr1 expr2
#'
#' @docType data
#'
#' @usage
#' data(expr1)
#' data(expr2)
#'
#' @format A matrix of integers, individuals as rows and genes as columns.
#'
#'
#' @keywords datasets
#'
#' @seealso \code{\link{genepos}}, \code{\link{f2cross}}, \code{\link{pmap}}
#'
#' @examples
#' data(expr1)
#' data(expr2)
#'
#' # identify the common individuals
#' id <- findCommonID(rownames(expr1), rownames(expr2))
#'
#' # correlation between tissues for each gene
#' rho <- corbetw2mat(expr1[id$first,], expr2[id$second,])
#' hist(rho, breaks=100)
"expr1"
