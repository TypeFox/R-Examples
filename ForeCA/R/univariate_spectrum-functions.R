#' @rdname mvspectrum
#' @keywords manip
#' @description
#' \code{get_spectrum_from_mvspectrum} extracts the spectrum of one time series from an
#' \code{"mvspectrum"} object by taking the i-th diagonal entry for each frequency.
#' @param which integer(s); the spectrum of which series whould be extracted. By default,
#' it returns all univariate spectra as a matrix (frequencies in rows).
#' @export
#' @return
#' \code{get_spectrum_from_mvspectrum} returns either a matrix of all univariate spectra,
#' or one single column (if \code{which} is specified.)
#' 
#' @examples
#' 
#' XX <- matrix(rnorm(1000), ncol = 2)
#' SS <- mvspectrum(XX, "direct")
#' ss1 <- mvspectrum(XX[, 1], "direct")
#' 
#' SS.1 <- get_spectrum_from_mvspectrum(SS, 1)
#' plot.default(ss1, SS.1)
#' abline(0, 1, lty = 2, col = "blue")
#' 

get_spectrum_from_mvspectrum <- function(mvspectrum.output, 
                                         which = seq_len(dim(mvspectrum.output)[2])) {
  num.series <- dim(mvspectrum.output)[2]
  
  stopifnot(all(which > 0),
            all(which <= num.series))
  
  if (is.null(num.series) || num.series == 1) {
    # just one dimensional series
    return(mvspectrum.output)
  } else {
    all.spectra <- t(apply(mvspectrum.output, 1, diag))
    tmp <- all.equal(rep(0, length(all.spectra)), c(Im(all.spectra)))
    if (!isTRUE(tmp)) {
      cat(tmp)
      warning("The multivariate spectrum has imaginary elements in the diagonal.", 
              " Please check your spectrum estimates again (and set 'inverse = FALSE'",
              " in 'normalize_mvspectrum' if you have used this function).")
    }
    return(Re(all.spectra[, which]))
  }
} 

#' @rdname mvspectrum
#' @keywords manip
#' @description
#' \code{spectrum_of_linear_combination} computes the spectrum of the linear
#' combination  \eqn{\mathbf{y}_t = \mathbf{X}_t \boldsymbol \beta} of \eqn{K} 
#' time series \eqn{\mathbf{X}_t}.  This can be efficiently computed by the 
#' quadratic form
#' \deqn{
#'   f_{y}(\lambda) = \boldsymbol \beta' f_{\mathbf{X}}(\lambda) \boldsymbol \beta \geq 0,
#' }
#' for each \eqn{\lambda}. This holds for any \eqn{\boldsymbol \beta} 
#' (even \eqn{\boldsymbol \beta = \boldsymbol 0} -- not only for 
#' \eqn{||\boldsymbol \beta ||_2 = 1}.
#' For \eqn{\boldsymbol \beta = \boldsymbol e_i} (the i-th basis vector) 
#' this is equivalent to \code{get_spectrum_from_mvspectrum(..., which = i)}.
#' 
#' @param beta numeric; vector \eqn{\boldsymbol \beta} that defines the linear
#' combination.
#' @export
#' @return
#' \code{spectrum_of_linear_combination} returns a vector with length equal to 
#' the number of rows of \code{mvspectrum.output}. 
#' @examples
#' 
#' XX <- matrix(arima.sim(n = 1000, list(ar = 0.9)), ncol = 4)
#' beta.tmp <- rbind(1, -1, 2, 0)
#' yy <- XX %*% beta.tmp
#' 
#' SS <- mvspectrum(XX, "wosa")
#' ss.yy.comb <- spectrum_of_linear_combination(SS, beta.tmp)
#' ss.yy <- mvspectrum(yy, "wosa")
#' 
#' plot(ss.yy, log = TRUE) # using plot.mvspectrum()
#' lines(ss.yy.comb, col = "red", lty = 1, lwd = 2) 
#' 
spectrum_of_linear_combination <- function(mvspectrum.output, beta) {
  
  stopifnot(dim(mvspectrum.output)[2] == length(beta))
  
  num.freqs <- dim(mvspectrum.output)[1]
  if (all(beta == 0)) {
    spec.dens.est <- rep(0, num.freqs)
  } else {
    spec.dens.est <- apply(mvspectrum.output, 1, quadratic_form, vec = beta)
    tmp <- all.equal(rep(0, length(spec.dens.est)), Im(spec.dens.est))
    
    if (!isTRUE(tmp)) {
      cat(tmp)
      warning("The linear combination of spectra has imaginary values.",
              " Please check your multivariate spectrum estimates again ")
    }
    spec.dens.est <- Re(spec.dens.est)
    
    # numerically sometimes this can be < 0; but this is just rounding error; set them to 0.
    # Adjust the overall vector so it has the same mean again as before setting it to 0.
    # Otherwise the average will be to large (as we remove negative values).
    ind.neg <- (spec.dens.est < 0)
    total.neg.values <- sum(spec.dens.est[ind.neg])
    spec.dens.est[ind.neg] <- 0
    spec.dens.est <- spec.dens.est * (1 + total.neg.values / sum(spec.dens.est))
  }
  return(spec.dens.est)
}