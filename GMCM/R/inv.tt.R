#' @rdname tt
#' @param par A vector of length 4 where \code{par[1]} is the mixture
#'   proportion, \code{tpar[2]} the mean, \code{tpar[3]} the standard deviation,
#'   and \code{tpar[4]} the correlation.
#' @return  \code{inv.tt} returns \code{tpar} as described above.
#' @keywords internal
inv.tt <- function(par, d, positive.rho) {
  if (par[3] < 0)
    stop("The standard deviation par[3] must be greater than zero.")
  if (par[1] < 0 | par[1] > 1)
    stop("The probability par[1] must in the interval (0,1).")
  if (par[2] < 0)
    stop("The mean par[2] must be greater than zero.")
  if ((par[4] < 0 | par[4] > 1) & positive.rho)
    stop("The correlation par[4] must be between 0 and 1.")
  if ((par[4] < -1/(d-1) | par[4] > 1) & !positive.rho)
    stop("The correlation par[4] must be between -1/(d-1) = ",
         round(-1/(d-1),3), " and 1.")
  tpar      <- NA
  tpar[1]   <- logit(par[1])
  tpar[2:3] <- log(par[2:3])
  if (positive.rho) {
    tpar[4] <- logit(par[4])
  } else {
    tpar[4] <- rho.transform(par[4], d)
  }
  return(tpar)
}
