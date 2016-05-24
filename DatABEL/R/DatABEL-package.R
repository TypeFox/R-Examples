#' DatABEL package for fast consecutive access to large out-of-RAM
#' stored matrices
#'
#' A package interfacing FILEVECTOR C++ library
#' for storage of and fast consecutive access to
#' large data matrices in out-of-RAM disk mode
#' with regulated cache size. Columns of matrix
#' are accessible very quickly.
#'
#' @seealso
#' \code{\link{apply2dfo}},
#' \code{\link{databel2matrix}},
#' \code{\link{databel2text}},
#' \code{\link{extract_text_file_columns}},
#' \code{\link{matrix2databel}},
#' \code{\link{text2databel}},
#' \code{\linkS4class{databel}}
#'
#' @docType package
#' @name DatABEL-package
#' @aliases "DatABEL"
#' @useDynLib DatABEL
#' @import methods
#'
#' @author Yurii Aulchenko (R code), Stepan Yakovenko (R and C++
#' code), Andrey Chernyh (C++ code)
#'
NULL
