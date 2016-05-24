## bssbsb-data.R

#' @name bssbsb
#' @title BSS/BSB backcross data
#'
#' @description Data from two densely genotyped backcrosses.
#'
#' @details
#' There are 94 individuals from each of two interspecific backcross: (C57BL/6J
#' \eqn{\times}{x} \emph{M. spretus}) \eqn{\times}{x} C57BL/6J and (C57BL/6J
#' \eqn{\times}{x} SPRET/Ei) \eqn{\times}{x} SPRET/Ei.  They were typed on 1372
#' and 4913 genetic markers, respectively, with 904 markers in common.
#'
#' These data are from September, 2000.  Updated data are available.
#'
#' @format An object of class \code{cross}.  See \code{\link[qtl]{read.cross}}
#' for details.
#' @references Rowe, L. B., Nadeau, J. H., Turner, R., Frankel, W. N., Letts,
#' V. A., Eppig, J. T., Ko, M. S., Thurston, S. J. and Birkenmeier, E. H.
#' (1994) Maps from two interspecific backcross DNA panels available as a
#' community genetic mapping resource.  \emph{Mamm. Genome} \bold{5}, 253--274.
#'
#' Broman, K. W., Rowe, L. B., Churchill, G. A. and Paigen, K. (2002) Crossover
#' interference in the mouse. \emph{Genetics} \bold{160}, 1123--1131.
#' @source Lucy Rowe, Jackson Laboratory
#' @keywords datasets
#' @examples
#'
#' data(bssbsb)
#' summary(bssbsb)
#' \dontrun{plot(bssbsb)}
#'
NULL
