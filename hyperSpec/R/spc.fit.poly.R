##' Polynomial Baseline Fitting
##' These functions fit polynomal baselines.
##' 
##' Both functions fit polynomials to be used as baselines. If \code{apply.to}
##' is \code{NULL}, a \code{hyperSpec} object with the polynomial coefficients
##' is returned, otherwise the polynomials are evaluated on the spectral range
##' of \code{apply.to}.
##' 
##' \code{spc.fit.poly} calculates the least squares fit of order
##' \code{poly.order} to the \emph{complete} spectra given in \code{fit.to}.
##' Thus \code{fit.to} needs to be cut appropriately.
##'
##' @rdname baselines
##' @concept baseline
##' @aliases spc.fit.poly spc.fit.poly.below
##' @param fit.to \code{hyperSpec} object on which the baselines are fitted
##' @param apply.to \code{hyperSpec} object on which the baselines are evaluted
##'   If \code{NULL}, a \code{hyperSpec} object containing the polynomial
##'   coefficients rather than evaluted baselines is returned.
##' @param poly.order order of the polynomial to be used
##' @return \code{hyperspec} object containing the baselines in the spectra
##'   matrix, either as polynomial coefficients or as polynomials evaluted on
##'   the spectral range of \code{apply.to}
##' @author C. Beleites
##' @seealso \code{vignette ("baseline", package = "hyperSpec")}
##' @keywords manip datagen
##' @export
##' @examples
##' 
##' \dontrun{vignette ("baseline", package = "hyperSpec")}
##' 
##' baselines <- spc.fit.poly(chondro[,, c (625 ~ 640, 1785 ~ 1800)], chondro)
##' plot(chondro - baselines, "spcprctl5")
##' 
spc.fit.poly <- function (fit.to, apply.to = NULL, poly.order = 1){
  chk.hy (fit.to)
  if (! is.null (apply.to))
    chk.hy (apply.to)

  validObject (fit.to)
  validObject (apply.to)

  x <- fit.to@wavelength
  x <- outer(x, 0 : poly.order, "^")             # Vandermonde matrix of x
  
  p <- apply (fit.to, 1, function (y, x){qr.solve (x, y)}, x)

  if (is.null (apply.to)){
    colnames (p@data$spc) <- paste ("x^", 0 : poly.order, sep="")

    p
  } else {
    wl <- apply.to@wavelength;
    x <- outer(wl, 0 : poly.order, "^")             # Vandermonde matrix of x
    apply.to@data$spc <- I (t (apply (p[[]], 1, function (p, x) {x %*% p}, x)))

    .wl(apply.to) <- wl
    colnames (apply.to@data$spc) <- format (wl, digits = 4)

    apply.to
  }
}

                                        
##' 
##' \code{spc.fit.poly.below} tries to fit the baseline on appropriate spectral
##' ranges of the spectra in \code{fit.to}.  For details, see the
##' \code{vignette ("baseline")}.
##' @rdname baselines
##' @param npts.min minmal number of points used for fitting the polynomial
##' @param noise noise level to be considered during the fit. It may be given
##'   as one value for all the spectra, or for each spectrum separately.
##' @export
##' @examples
##' 
##' baselines <- spc.fit.poly.below(chondro)
##' plot(chondro - baselines, "spcprctl5")
##' 
spc.fit.poly.below <- function (fit.to, apply.to = fit.to, poly.order = 1,
                                npts.min = NULL, noise = 0){
  chk.hy (fit.to)
  if (! is.null (apply.to))
    chk.hy (apply.to)

  validObject (fit.to)
  validObject (apply.to)

  if (is.null (npts.min)){
    npts.min <- max (round (nwl(fit.to) * 0.05), 3 * (poly.order + 1))
    cat ("Fitting with npts.min = ",  npts.min, "\n")
  } else  if (npts.min <= poly.order){
    npts.min <- poly.order + 1
    warning (paste ("npts.min too small: adjusted to", npts.min))
  }

  if (length (noise) == 1)
    noise <- rep (noise, nrow (fit.to))

  vdm <- outer(fit.to@wavelength, 0 : poly.order, "^")
  y <- t(fit.to [[]])

  p <- matrix (nrow = nrow(fit.to) , ncol = poly.order + 1)
  for (i in row.seq (fit.to)){
    use.old <- logical (nwl (fit.to))
    use <- !use.old

    repeat {
      p[i,] <- qr.solve (vdm[use,], y[use, i])
      bl <- vdm %*% p [i,]
      use.old <- use
      use <- y[, i] < bl + noise [i]
      if ((sum (use, na.rm=TRUE) < npts.min) || all (use == use.old, na.rm = TRUE))
        break
    }
  }
  if (is.null (apply.to)){
    fit.to@data$spc <- p
    .wl (fit.to) <- 0 : poly.order
    colnames (fit.to@data$spc) <- paste ("x^", 0 : poly.order, sep="")

    validObject (fit.to)
    
    fit.to
  } else {
    x <- apply.to@wavelength

    vdm <- outer(x, 0 : poly.order, "^")             # Vandermonde matrix of x

    apply.to@data$spc <- I (t (apply (p, 1, function (p, x) {x %*% p}, vdm)))

    .wl(apply.to) <- x
    colnames (apply.to@data$spc) <- format (x, digits = 4)

    validObject (apply.to)
    
    apply.to
  }
}
