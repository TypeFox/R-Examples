#' @name gvc
#' @docType package
#' @title Global Value Chain analysis
#' @author
#' Bastiaan Quast \email{bquast@@gmail.com}
#' Victor Kummritz
#' @description Several tools for Global Value Chain ('GVC') analysis are implemented.
#' @seealso http://qua.st/decompr
#' @references {Wang, Zhi, Shang-Jin Wei, and Kunfu Zhu. Quantifying international production sharing at the bilateral and sector levels. No. w19677. National Bureau of Economic Research, 2013.}
NULL
.onAttach <- function(...) {
  packageStartupMessage('If you use gvc for data analysis,
please cite both R and gvc,
using citation() and citation("gvc") respectively.
                        ')}
