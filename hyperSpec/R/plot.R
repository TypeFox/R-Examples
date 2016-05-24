###-----------------------------------------------------------------------------
###
###  plot methods
###

###-----------------------------------------------------------------------------
###
### .plot: main switchyard for plotting functions
###

.plot <-  function (x, y, ...){
  ##    'spc'        ... spectra
  ##    'map'        ... map
  ##    'voronoi'    ... voronoi tiled map
  ##    'mat'        ... spectra matrix
  ##    'c'          ... concentration: plotc
  ##    'ts'         ... time series: plotc
  ##    'depth'      ... concentration or time series
  ##    'spcmeansd'  ... mean spectrum +- 1 standard deviation
  ##    'spcprctile' ... median spectrum , 16th and 84th percentile
  ##    'spcprctl5'  ... spcprctile plus 5th and 95th percentile
  
  dots <- list(...)          # to allow optional argument checks
  
  if (missing (y)){
    stop ("second argument to plot is missing. Should be a character indicating the type of plot.")
    y = "spc"
  }
  
  switch (tolower (y),
          spc = plotspc (x, ...),
          spcmeansd = {
            dots <- modifyList (list (object = mean_pm_sd (x),
                                      fill = c (1, NA, 1)
                                      ),
                                dots)
            do.call (plotspc, dots)
          },
          spcprctile = {
            dots <- modifyList (list (object = quantile (x, probs = c (0.16, 0.5, 0.84)),
                                      fill = c (1, NA, 1)
                                      ),
                                dots)
            do.call (plotspc, dots)
          },
          spcprctl5 = {
            dots <- modifyList (list (object = quantile (x, probs = c (0.05, 0.16, 0.5, 0.84, 0.95)),
                                      fill = c (1, 2, 3, 2, 1),
                                      fill.col = c("#00000040")
                                      ),
                                dots)
            do.call (plotspc, dots)
          },
          map = plotmap (x, ...),
          voronoi = plotvoronoi (x, ...),
          mat = plotmat (x, ...),
          c = plotc (x, ...),
          ts = plotc (x, spc ~ t, ...),
          depth = plotc (x, spc ~ z, ...),
          stop (paste ("y = ", y, "unknown.", collapse = " "))
          )
}

##' @noRd
##' @export
setGeneric ('plot')

##' Plotting hyperSpec Objects
##' 
##' Plotting \code{hyperSpec} objects. The \code{plot} method for
##' \code{hyperSpec} objects is a switchyard to \code{\link{plotspc}},
##' \code{\link{plotmap}}, and \code{\link{plotc}}.
##' 
##' It also supplies some convenient abbrevations for much used plots.
##' 
##' If \code{y} is missing, \code{plot} behaves like \code{plot (x, y =
##' "spc")}.
##' 
##' Supported values for \code{y} are:
##' 
##' \describe{ \item{"spc"}{calls \code{\link{plotspc}} to produce a spectra
##' plot.}
##' 
##' \item{"spcmeansd"}{plots mean spectrum +/- one standard deviation}
##' 
##' \item{"spcprctile"}{plots 16th, 50th, and 84th percentile spectre. If the
##' distributions of the intensities at all wavelengths were normal, this would
##' correspond to \code{"spcmeansd"}. However, this is frequently not the case.
##' Then \code{"spcprctile"} gives a better impression of the spectral data
##' set.}
##' 
##' \item{"spcprctl5"}{like \code{"spcprctile"}, but additionally the 5th and
##' 95th percentile spectra are plotted.}
##' 
##' \item{"map"}{calls \code{\link{plotmap}} to produce a map plot.}
##' 
##' \item{"voronoi"}{calls \code{\link{plotvoronoi}} to produce a Voronoi plot
##' (tesselated plot, like "map" for hyperSpec objects with uneven/non-rectangular
##' grid).}
##'
##' \item{"mat"}{calls \code{\link{plotmat}} to produce a plot of the spectra
##' matrix (not to be confused with \code{\link[graphics]{matplot}}).}
##' 
##' \item{"c"}{calls \code{\link{plotc}} to produce a calibration (or time
##' series, depth-profile, or the like)}
##' 
##' \item{"ts"}{plots a time series: abbrevation for \code{\link{plotc} (x,
##' use.c = "t")}}
##' 
##' \item{"depth"}{plots a depth profile: abbrevation for \code{\link{plotc}
##' (x, use.c = "z")}} }
##' 
##' @name plot-methods
##' @rdname plot
##' @aliases plot plot,ANY,ANY-method plot,hyperSpec,character-method
##'   plot,hyperSpec,missing-method
##' @docType methods
##' @param x the \code{hyperSpec} object
##' @param y selects what plot should be produced
##' @param ... arguments passed to the respective plot function
##' @author C. Beleites
##' @seealso \code{\link{plotspc}} for spectra plots (intensity over
##'   wavelength),
##' 
##' \code{\link{plotmap}} for plotting maps, i.e. color coded summary value on
##'   two (usually spatial) dimensions.
##' 
##' \code{\link{plotc}}
##' 
##' \code{\link[graphics]{plot}}
##' @keywords methods hplot
##' @export
##' @examples
##' 
##' plot (flu)
##' 
##' plot (flu, "c")
##' 
##' plot (laser, "ts")
##' 
##' spc <- apply (chondro, 2, quantile, probs = 0.05)
##' spc <- sweep (chondro, 2, spc, "-")
##' plot (spc, "spcprctl5")
##' plot (spc, "spcprctile")
##' plot (spc, "spcmeansd")
##'
### use plotspc as default plot function
setMethod ("plot",
           signature (x = "hyperSpec", y = "missing"),
           function (x, y, ...) plotspc (x, ...)
           )

### allow choice of spectral or map plot by second argument
##' @rdname plot
##' @export
setMethod ("plot",
           signature (x = "hyperSpec", y = "character"), .plot)
