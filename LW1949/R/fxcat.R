#' Define Effect Category
#'
#' Define three effect categories, 0 for none affected, 100 for all affected,
#'   and 50 for other proportions affected.
#' @param dat
#'   A data frame of toxicity data, including at least two variables:
#'     ntot (the number of individuals tested) and
#'     nfx (the number of affected individuals).
#' @return
#'   An integer vector the same length as \code{prob} with
#'     categories of 0, 50, or 100.
#' @export
#' @examples
#' toxdat <- data.frame(
#'   dose=c(0.0625, 0.125, 0.25, 0.5),
#'   ntot=rep(8, 4),
#'   nfx = c(0, 4, 6, 8))
#' cbind(toxdat, fxcat(toxdat))

fxcat <- function(dat) {
  if (!is.data.frame(dat)) stop("Input must be a data frame.")
  if (any(is.na(match(c("ntot", "nfx"), names(dat))))) {
    stop("Input must include at least two variables: ntot, nfx")
  }
  categ <- rep(50, dim(dat)[1])
  categ[with(dat, !is.na(nfx) & nfx==0)] <- 0
  categ[with(dat, !is.na(ntot) & !is.na(nfx) & ntot==nfx)] <- 100
  categ[with(dat, is.na(ntot) | is.na(nfx))] <- NA
  categ[with(dat, (!is.na(nfx) & !is.na(ntot) & nfx>ntot) |
      (!is.na(nfx) & nfx<0) |
      (!is.na(ntot) & ntot<1))] <- NA
  as.integer(categ)
  }
