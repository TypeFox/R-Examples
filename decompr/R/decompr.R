#' @name decompr
#' @docType package
#' @title Export Decomposition using the Wang-Wei-Zhu and Leontief decompositions algorithms.
#' @author
#' Bastiaan Quast \email{bquast@@gmail.com}
#' Fei Wang
#' Victor Kummritz
#' @seealso http://qua.st/decompr
#' @references {Wang, Zhi, Shang-Jin Wei, and Kunfu Zhu. Quantifying international production sharing at the bilateral and sector levels. No. w19677. National Bureau of Economic Research, 2013.}
NULL
#' @name final
#' @docType data
#' @title Leather Example
#' @description the final demand data
NULL
#' @name inter
#' @docType data
#' @title Leather Example
#' @description the intermediate demand data
NULL
#' @name countries
#' @docType data
#' @title Leather Example
#' @description the names of the countries data
NULL
#' @name industries
#' @docType data
#' @title Leather Example
#' @description the names of the industries data
NULL
#' @name out
#' @docType data
#' @title Leather Example
#' @description final output
NULL
.onAttach <- function(...) {
  packageStartupMessage('If you use decompr for data analysis,
please cite both R and decompr,
using citation() and citation("decompr") respectively.
')}
