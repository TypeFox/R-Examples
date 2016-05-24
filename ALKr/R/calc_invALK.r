#' Classic ALK
#' 
#' Returns an inverse Age-Length Key calculated from a matrix with the count of
#' individuals per age- and length-class, as described by Clark (1981), Bartoo
#' and Parker (1983) and Hilborn and Walters (1992)
#' 
#' @param x A matrix with \code{i} lines and \code{j} columns, where
#' \code{x[i, j]} is the count of individuals of length \code{i} and age
#' \code{j}.
#' @param fi A vector of length \code{i} where \code{fi[i]} is the number of
#' fish in the length-class \code{j} on the population from which \code{x} was
#' sampled.
#' 
#' @return A matrix with the probability of an individual of age \code{j}
#' having length \code{i}, i.e. \eqn{P(i|j)}.
#' 
#' @references
#' Bartoo, N.W., Parker, K.R. (1983). Stochastic age-frequency estimation using
#' the von Bertalanffy growth equation. \emph{Fishery Bulletin}, \strong{81}/1,
#' 91-96
#' 
#' Clark, W.G. (1981). Restricted Least-Squares Estimates of Age Composition
#' from Length Composition. \emph{Canadian Journal of Fisheries and Aquatic
#' Sciences}, \strong{38}/3, 297-307. DOI: \code{10.1139/f81-041}
#' 
#' Hilborn, R., Walters, C.J. (1992). Quantitative fisheries stock assessment:
#' Choice, dynamics and uncertainty. \emph{Reviews in Fish Biology and
#' Fisheries}, \strong{2}/2, 177-178. DOI: \code{10.1007/BF00042883}
#' 
#' @examples
#' data(hom)
#' calc_invALK(hom$otoliths[[1]], hom$F1992)
#' 
#' @export
calc_invALK <- function(x, fi) {
  N <- fi * calc_ALK(x)
  denom <- colSums(N)
  denom[denom==0] <- 1
  t(t(N) / denom)
}

