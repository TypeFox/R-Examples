#
# All code copyright 2013 Timothy H. Keitt
#

#' Little Rock Lake zooplankton dataset
#' 
#' Contains time series for 10 dominant crustaceous species of zooplanking
#' sampled from Little Rock Lake, Wisconsin. Samples come from two basins: one
#' treated to lower pH and the other an untreated reference.
#'
#' @name lrlake
#' @docType data
#' @format A data frame with 592 observations on the following 18 variables.
#' @references Keitt, T. H. 2008. Coherent ecological dynamics induced by
#' large-scale disturbance. Nature 454:331-4. doi:10.1038/nature06935.
#' @source http://lter.limnology.wisc.edu/
#' @keywords datasets
#' @examples
#' 
#' data(lrlake)
#' x = subset(lrlake, Basin == "Reference", LRL.Day)
#' y = subset(lrlake, Basin == "Reference", -(1:8))
#' matplot(x, y, type = "l", lty = 1)
#' 
NULL





#' Wavelet transform of multivariate time series
#' 
#' Computes continuous wavelet transform of multiple irregularly sampled time
#' series.
#' 
#' \tabular{ll}{ Package: \tab mvcwt\cr Type: \tab Package\cr Version: \tab
#' 1.3\cr Date: \tab 2013-10-27\cr License: \tab GPL\cr } The main functions
#' are \code{\link{mvcwt}}, which computes the wavelet transform of multiple
#' time series, and \code{\link{wmr}}, which computes the wavelet modulus
#' ratio, a measure of time series coherence.
#' 
#' Note that this is a complete rewrite of the code used in the reference
#' below, and as such it is not well tested. It may give different or
#' inaccurate results. I recommend you run tests on known data.
#' 
#' The most recent development version of this package can be found at
#' \url{https://bitbucket.org/tkeitt/mvcwt/overview}.
#' 
#' @name mvcwt-package
#' @docType package
#' @author Timothy H. Keitt (\url{http://www.keittlab.org})
#' 
#' Tim Keitt <tkeitt@@gmail.com>
#' @references Keitt, T. H. 2008. Coherent ecological dynamics induced by
#' large-scale disturbance. Nature 454:331-4. doi:10.1038/nature06935.
#' @keywords package
#' @examples
#' 
#' \dontrun{
#'   x = seq(-pi, pi, len = 200)
#'   y1 = sin(8 * x) + sin(32 * x)
#'   y2 = sin(8 * (x + pi/8)) + sin(32 * x)
#'   matplot(x, cbind(y1, y2), type = "l", lty = 1)
#'   w = mvcwt(x, cbind(y1, y2))
#'   plot(w, var = 1:2, scale = 2^seq(log2(min(w$y)), log2(max(w$y)), len = 5))
#'   mr = wmr(w, smoothing = 2)
#'   image(mr, reset.par = FALSE)
#'   contour(mr, levels = c(0.01, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.99), add = TRUE)}
#'   
NULL

#' @import foreach
#' @import RColorBrewer
#' @import grDevices
NULL

#' Utility functions
#' 
#' Functions for calculating scales and locations to analyze
#' 
#' @param x a vector of values
#' @param min smallest value in sequence
#' @param max largest value in sequence
#' @param nbins how many intervals
#' @param nsteps length of returned sequence
#' 
#' @return
#' \code{get.nscales}: length of \code{x} \code{\link{unlist}}ed and \code{\link{as.vector}}ized.
#' 
#' \code{get.min.scale}: twice the median distance between successive values of \code{x}
#' 
#' \code{get.max.scale}: 1/2 the maximum distance between values in \code{x}
#' 
#' \code{log2Bins}: a sequence of values on a log2 scale
#' 
#' \code{regularize}: a regular sequence of values
#' 
#' @rdname utils
#' @export
get.nscales = function(x) length(as.vector(unlist(x)))

#' @rdname utils
#' @export
get.min.scale = function(x) 2 * median(diff(as.vector(unlist(x))))

#' @rdname utils
#' @export
get.max.scale = function(x) diff(range(as.vector(unlist(x)))) / 2

#' @rdname utils
#' @export
log2Bins = function(min, max, nbins)
{
  2 ^ midp(seq(log2(min), log2(max), length = nbins + 1))
}

#' @rdname utils
#' @export
regularize = function(x, nsteps = length(as.vector(unlist(x))))
{
  seq(min(x), max(x), length = nsteps)
}

#' Computes the wavelet transform of a multivariate time series
#' 
#' This function takes set a set of seqences as columns in a matrix and
#' computes the continuous wavelet transform on each.
#' 
#' @aliases mvcwt print.mvcwt
#' @param x sample locations
#' @param y one or more columns of samples corresponding to \code{x}
#' @param scale.exp scale output
#' @param nscales number of scales to analyze
#' @param min.scale minimum scale in units of \code{x}
#' @param max.scale maximum scale in units of \code{x}
#' @param scales a set of scales to analyze; overrides all other scale arguments
#' @param loc the loci at which to evalues the wavelet function
#' @param wave.fun a wavelet function
#' 
#' @export mvcwt
mvcwt = function(x, y,
                 scale.exp = 0.5,
                 nscales = get.nscales(x),
                 min.scale = get.min.scale(x),
                 max.scale = get.max.scale(x),
                 scales = log2Bins(min.scale, max.scale, nscales),
                 loc = regularize(x), wave.fun = "Morlet")
{
  s = 1 # this is a workaround for a bug in R's code checking
  wave.fun = match.fun(wave.fun)
  x = as.vector(unlist(x))
  lmat = lagMat(x, loc)
  y = matrix(unlist(y), nrow = length(x))
  w = foreach(s = scales) %dopar%
  {
    Conj(wave.fun(lmat / s) / s ^ scale.exp) %*% y
  }
  w = array(unlist(w), dim = c(length(loc), ncol(y), length(scales)))
  w = aperm(w, perm = c(1, 3, 2))
  structure(list(x = loc, y = scales, z = w), class = "mvcwt")
}

#' @export
print.mvcwt = function(x, ...)
{
  print(str(x), ...)
  invisible(x)
}

#' Compute the wavelet modulus ratio of multivariate data
#' 
#' Computes the wavelet modulus ratio described in Keitt (2008). A value of one
#' indicated perfect synchrony among time series and a value of zero, perfect
#' compensation.
#' 
#' @param w an object such as returned by \code{\link{mvcwt}}
#' @param smoothing width of smoothing kernel; larger values give more
#' smoothing
#' 
#' @return an object of class "mvcwt"
#' 
#' @author Timothy H. Keitt
#' 
#' @seealso \code{\link{mvcwt}}, \code{\link{image.mvcwt}}
#' 
#' @references Keitt, T. H. 2008. Coherent ecological dynamics induced by
#' large-scale disturbance. Nature 454:331-4. doi:10.1038/nature06935.
#' 
#' @keywords wavelets
#' 
#' @examples
#' \dontrun{
#' data(lrlake)
#' x = subset(lrlake, Basin == "Treatment", LRL.Day) / 365.25
#' y = subset(lrlake, Basin == "Treatment", -(1:8))
#' w = mvcwt(x, y, min.scale = 0.25, max.scale = 4)
#' mr = wmr(w)
#' image(mr, reset.par = FALSE)
#' contour(mr, bound = NA, add = TRUE)}
#' 
#' @export wmr
wmr = function(w, smoothing = 1)
{
  with(w, {
    mods = Mod(rowSums(z, dims = 2))
    smod = rowSums(Mod(z), dims = 2)
    lmat = lagMat(x)
    exports = c("Gauss", "smoothing")
    flibs = c("mvcwt")
    modrat = foreach(i = 1:length(y),
                     .combine = cbind,
                     .export = exports,
                     .packages = flibs) %dopar%
    {
      kern = Gauss(lmat / y[i] / smoothing)
      modsv = kern %*% mods[,i]
      smodv = kern %*% smod[,i]
      modsv / smodv
    }
    dim(modrat) = c(length(x), length(y), 1)
    structure(list(x = x, y = y, z = modrat), class = "mvcwt")
  })
}

#' Boot strap p-values for wavelet modulus ratio
#' 
#' Performs a phase-randomization bootstrap estimate of the null hypothesis of
#' independent time series
#' 
#' The phases are randomized \code{reps} times for each combination of input
#' variable and scale. This package depends heavily on the \code{dopar}
#' function in the \code{foreach} package. If you do not have a lot of cores
#' available to you, you may need to let this run overnight.
#' 
#' @param w an object such as returned by \code{\link{mvcwt}}
#' @param smoothing degree of smoothing; larger values give greater smoothing
#' @param reps number of repetitions
#' @param mr.func a function taking a "mvcwt" object to be applied to each
#' trial
#' 
#' @return an object of class "mvcwt" suitable for use with
#' \code{\link{contour.mvcwt}}.
#' 
#' @author Timothy H. Keitt
#' 
#' @seealso \code{\link{mvcwt}}, \code{\link{wmr}}
#' 
#' @references Keitt, T. H. 2008. Coherent ecological dynamics induced by
#' large-scale disturbance. Nature 454:331-4. doi:10.1038/nature06935.
#' 
#' @keywords statistics
#' @export wmr.boot
wmr.boot = function(w, smoothing = 1, reps = 1000, mr.func = "wmr")
{
  mr.func = match.fun(mr.func)
  mr.obs = mr.func(w, smoothing = smoothing)
  with(w, {
    nloc = length(x)
    nvars = dim(z)[3]
    nscales = length(y)
    exports = c("reps", "wmr", "lagMat", "regularize", "mr.func",
                "Gauss", "mr.obs", "nscales", "smoothing")
    flibs = c("mvcwt")
    mr.obs$z.boot = foreach(i = 1:nscales,
                            .combine = c,
                            .export = exports,
                            .packages = flibs) %dopar%
    {
      mr.boot = foreach(j = 1:reps,
                        .combine = cbind,
                        .inorder = FALSE) %dopar%
      {
        rphase = t(array(runif(nvars, -pi, pi), dim = c(nvars, nloc)))
        zp = z[, i,, drop = FALSE] * complex(argument = rphase)
        as.vector(mr.func(list(x = x, y = y[i], z = zp), smoothing = smoothing)$z)
      }
      res = foreach(j = 1:nloc, .combine = c) %dopar%
      {
         ecdf(mr.boot[j, ])(mr.obs$z[j, i,])
      }
      res
    }
    dim(mr.obs$z.boot) = c(length(x), length(y), 1)
    return(mr.obs)
  })
}

#' Plot wavelet output
#' 
#' Plot multivariate wavelet output across variables, scales or both.
#' 
#' Makes one or more plots on the graphics device. Total number of plots is
#' limited to 10.
#' 
#' @param x an object such as produced by \code{\link{mvcwt}}
#' @param var which variables to plot; can be a vector
#' @param scale which scales to plot; can be a vector; closest scale is picked
#' @param titles plot titles on each sub-plot?
#' @param z.fun apply function to data prior to plotting
#' @param \dots additional graphical parameters passed to \code{\link{plot}}
#' 
#' @return \code{x} is returned invisibly
#' 
#' @author Timothy H. Keitt
#' 
#' @seealso \code{\link{mvcwt}}
#' 
#' @keywords graphics
#' 
#' @examples
#' \dontrun{
#' data(lrlake)
#' x = subset(lrlake, Basin == "Reference", LRL.Day) / 365.25
#' y = subset(lrlake, Basin == "Reference", -(1:8))
#' w = mvcwt(x, y, min.scale = 0.25, max.scale = 4)
#' plot(w, var = 1:10)}
#' 
#' @export plot.mvcwt
plot.mvcwt = function(x, var = 1, scale = 1, titles = TRUE, z.fun = "Re", ...)
{
  z.fun = match.fun(z.fun)
  opar = par(no.readonly = TRUE)
  on.exit(par(opar))
  with(x, {
    nvar = length(var)
    nscal = length(scale)
    par(mfrow = c(nvar * nscal, 1), mar = rep(0.6, 4), oma = rep(5, 4), xpd = NA)
    for ( s in scale )
    {
      scale.i = which.min(abs(y - s))
      for ( i in 1:nvar )
      {
        z.out = z.fun(z[, scale.i, var[i]])
        plot(x, z.out, xlab = NA, ylab = NA, axes = FALSE, type = "l", lwd = 2, ...)
        lines(range(x), rep(median(z.out), 2), lty = 2)
        if (titles )
        {
          tstr = paste("Scale =", signif(y[scale.i], 2), "Var =", var[i], sep = " ")
          title(main = tstr, line = 0.25)
        }
        if ( i %% 2 ) axis(2) else axis(4)
      }
    }
    axis(1)
    mtext("Location", 1, 3, outer = TRUE)
    mtext("Value", 2, 3, outer = TRUE)
  })  
  invisible(x)
}

#' Draw a heatmap of a \code{\link{mvcwt}} object
#' 
#' Draws one or more heatmaps
#' 
#' This function will draw a series of heatmaps on the graphical device. If you
#' want to add additional graphical elements, set \code{reset.par} to false.
#' 
#' @param x an object as returned by \code{\link{mvcwt}}
#' @param z.fun a function applied to the data before plotting
#' @param bound if finite, draw lines \code{bound * scale} units inside the
#' plot boundaries
#' @param reset.par if true, reset graphical parameters on exit
#' @param \dots additional arguments passed to \code{\link{image}}
#' 
#' @return \code{x} is returned invisibly
#' 
#' @author Timothy H. Keitt
#' 
#' @seealso \code{\link{mvcwt}}, \code{\link{wmr}}
#' 
#' @keywords graphics
#' 
#' @examples
#' \dontrun{
#' data(lrlake)
#' x = subset(lrlake, Basin == "Treatment", LRL.Day) / 365.25
#' y = subset(lrlake, Basin == "Treatment", -(1:8))
#' w = mvcwt(x, y, min.scale = 0.25, max.scale = 4)
#' image(w, z.fun = "Mod")}
#' 
#' @export image.mvcwt
image.mvcwt = function(x, z.fun = "Re", bound = 1, reset.par = TRUE, ...)
{
  z.fun = match.fun(z.fun)
  opar = par(no.readonly = TRUE)
  if ( reset.par ) on.exit(par(opar))
  pal = colorRampPalette(rev(brewer.pal(11, 'Spectral')))(1024)
  with(x, {
    nvar = ifelse(length(dim(z)) == 3, dim(z)[3], 1)
    par(mfrow = c(nvar, 1), mar = rep(0.2, 4), oma = rep(5, 4))
    for ( i in 1:nvar )
    {
      image(x, y, z.fun(z[,,i]), log = "y", col = pal, axes = FALSE, ...)
      if ( i %% 2 ) axis(2) else axis(4)
      if ( exists("z.boot") && !is.null(z.boot) )
      {
        z.boot = 1 - abs(1 - 2 * z.boot)
        contour(x, y, z.boot[,,i], levels = 0.05, lty = 3, add = TRUE, drawlabels = FALSE)
        zb = p.adjust(as.vector(z.boot), method = "BY")
        dim(zb) = dim(z.boot)
        contour(x, y, zb[,,i], levels = 0.05, lwd = 2, add = TRUE, drawlabels = FALSE)
      }
      if ( is.finite(bound) )
      {
        lines(min(x) + bound * y, y, lty = 2, lwd = 2, col = "darkgrey")
        lines(max(x) - bound * y, y, lty = 2, lwd = 2, col = "darkgrey")
      }
      box()
    }
    axis(1)
    mtext("Location", 1, 3, outer = TRUE)
    mtext("Scale", 2, 3, outer = TRUE)
  })
  return(invisible(x))     
}

#' Make contour plot of a \code{\link{mvcwt}} object
#' 
#' Draws a contour plot
#' 
#' Draws a contour plot. If you want to add more plotting elements, set
#' \code{reset.pars} to false.
#' 
#' @param x an object produced by \code{\link{mvcwt}} or \code{\link{wmr}}
#' @param z.fun a function to apply to the data prior to plotting
#' @param bound if finite, draw boundary lines \code{bound * scale} from plot
#' boundaries
#' @param reset.par if true, reset the graphical parameters on exit
#' @param \dots passed to the \code{\link{contour}} function
#' 
#' @return The object \code{x} is returned invisibly.
#' 
#' @author Timothy H. Keitt
#' 
#' @seealso \code{\link{mvcwt}}, \code{\link{wmr}}
#' 
#' @keywords graphics
#' 
#' @examples
#' \dontrun{
#' data(lrlake)
#' x = subset(lrlake, Basin == "Treatment", LRL.Day) / 365.25
#' y = subset(lrlake, Basin == "Treatment", -(1:8))
#' w = mvcwt(x, y, min.scale = 0.25, max.scale = 4)
#' mr = wmr(w)
#' contour(mr)}
#'   
#' @export contour.mvcwt
contour.mvcwt = function(x, z.fun = "Re", bound = 1, reset.par = TRUE, ...)
{
  z.fun = match.fun(z.fun)
  opar = par(no.readonly = TRUE)
  if ( reset.par ) on.exit(par(opar))
  with(x, {
    nvar = ifelse(length(dim(z)) == 3, dim(z)[3], 1)
    par(mfrow = c(nvar, 1), mar = rep(0.2, 4), oma = rep(5, 4))
    for ( i in 1:nvar )
    {
      contour(x, y, z.fun(z[,,i]), log = "y", axes = FALSE, ...)
      if ( i %% 2 ) axis(2) else axis(4)
      if ( is.finite(bound) )
      {
        lines(min(x) + bound * y, y, lty = 2, lwd = 2, col = "darkgrey")
        lines(max(x) - bound * y, y, lty = 2, lwd = 2, col = "darkgrey")
      }
      box()
    }
    axis(1)
    mtext("Location", 1, 3, outer = TRUE)
    mtext("Scale", 2, 3, outer = TRUE)
  })
  return(invisible(x))     
}

midp <- function(x)
{
  (x[-1] + x[-length(x)]) / 2
}

lagMat <- function(from, to = regularize(from))
{
  outer(to, from, "-")
}

Gauss <- function(lag)
{
  exp(-lag ^ 2 / 2) / sqrt(2 * pi)
}
  
#' The Morlet function
#' 
#' Given a sequence of lag distances, this function returns the Morlet wavelet
#' 
#' This version of the Morlet is scaled so that the central frequency is
#' exactly 2Pi radians. This is the simple version of the Morlet, sometimes
#' referred to as a psuedo-wavelet as it it not precisely normalized, leading
#' to some leakage into the DC component. It is therefore unsuited to
#' reconstruction using the inverse transform.
#' 
#' @param lag A sequence of lag distances, typically a matrix
#' @return A set of Morlet filters, typically as a matrix
#' @author Timothy H. Keitt
#' @seealso \code{\link{mvcwt}}
#' @keywords wavelets
#' @examples
#' x = seq(-pi, pi, len = 256)
#' plot(x, Re(Morlet(x)), col = "darkblue", type = "l")
#' lines(x, Im(Morlet(x)), col = "darkred")
#' lines(range(x), rep(0, 2), lty = 2)
#' 
#' @export Morlet
Morlet <- function(lag)
{
  exp(-lag ^ 2 / 2 + 2i * pi * lag) / pi ^ 4
}

