##' Identifying Spectra and Spectral Data Points
##' This function allows to identify the spectrum and the wavelength of a point
##' in a plot produced by \code{\link{plotspc}}.
##' 
##' This function first finds the spectrum with a point closest to the clicked
##' position (see \code{\link[graphics]{locator}}). The distance to the clicked
##' point is evaluated relative to the size of the tolerance window.
##' 
##' In a second step, \code{max.fn} searches for the actual point to label
##' within the specified wavelength window of that spectrum. This allows to
##' label maxima (or minima) without demanding too precise clicks. Currently,
##' the following functions to determine the precise point: \tabular{ll}{
##' spc.point.default \tab uses the clicked wavelength together with its
##' spectral intensity\cr spc.point.max \tab the point with the highest
##' intensity in the wavelength window \cr spc.point.min \tab the point with
##' the lowest intensity in the wavelength window \cr spc.point.sqr \tab
##' maximum of a parabola fit throug the point with highest intensity and the
##' two surrounding points \cr } \code{point.fn} is called with the arguments
##' \code{wl} containing the considered wavelength window, \code{spc} the
##' respective intensities of the closest spectrum, and \code{wlclick} the
##' wavelength that was clicked. They return a vector of two elements
##' (wavelength and intensity).
##' 
##' As a last step, a label for the point produced by \code{formatter} and
##' plotted using \code{\link[graphics]{text}}. Currently, the following
##' \code{formatter}s are available: \tabular{ll}{ spc.label.default \tab
##' spectrum number, wavelength \cr spc.label.wlonly \tab wavelength\cr }
##' \code{formatter} functions receive the number of the spectrum \code{ispc},
##' the wavelength \code{wl}, and the spectral intensity \code{spc} and produce
##' a character variable suitable for labelling. The predefined formatters
##' surround the label text by spaces in order to easily have an appropriate
##' offset from the point of the spectrum.
##' 
##' The warning issued if no spectral point is inside the tolerance window may
##' be switched of by \code{warn = FALSE}. In that case, the click will produce
##' a row of \code{NA}s in the resulting data.frame.
##' 
##' \code{spc.identify} uses option \code{debuglevel} to determine whether debugging output should be
##' produced. \code{debuglevel == 2} will plot the tolerance window for every clicked point,
##' \code{debuglevel == 1} will plot the tolerance window only if no data point was inside.  See
##' \code{\link[hyperSpec:options]{hyperSpec options}} for details about retrieving and setting
##' options.
##' 
##' You may want to adjust the plot's \code{ylim} to ensure that the labels are
##' not clipped. As a dirty shortcut, \code{xpd = NA} may help.
##' 
##' @aliases spc.identify spc.label.default spc.label.wlonly spc.point.default
##'   spc.point.max spc.point.min spc.point.sqr
##' @param x either the abscissa coordinates or the list returned by
##'   \code{\link{plotspc}}
##' @param y the ordinate values. Giving \code{y} will override any values from
##'   \code{x$y}.
##' @param wavelengths the wavelengths for the data points. Giving
##'   \code{wavelengths} will override any values from \code{x$wavelengths}.
##' @param tol.wl,tol.spc tolerance in wavelength and spectral intensity to
##'   search around the clicked point. See details.
##' @param point.fn \code{function (wl, spc, wlclick)} to determine the actual
##'   point to label, see details.
##' @param formatter \code{function (i, wl, spc)} that produces the labels. If
##'   \code{NULL}, no labels are displayed.
##' @param ... passed to \code{\link[graphics]{text}} in order to produce the
##'   labels
##' @param  cex,adj,srt see \code{\link[graphics]{par}}
##' @param warn Should the user be warned if no point is in the considered
##'   window? In addition, see the discussion of option \code{debuglevel} in
##'   the details.
##' 
##' If \code{FALSE}, the resulting data.frame will have a row of \code{NA}s
##'   instead.
##' @param delta \code{spc.point.sqr} fits the parabola in the window wlclick
##'   \eqn{\pm}{+-} delta points.
##' @return a \code{data.frame} with columns \item{ispc}{spectra indices of the
##'   identified points, i.e. the rows of the \code{hyperSpec} object that was
##'   plotted.
##' 
##' If \code{ispc} is given, \code{ispc [i]} is returned rather than \code{i}.
##'   } \item{wavelengths}{the wavelengths of the identified points}
##'   \item{spc}{the intensities of the identified points}
##' @author C. Beleites
##' @seealso \code{\link[graphics]{locator}}, \code{\link{plotspc}},
##'   \code{\link[hyperSpec:options]{hyperSpec options}}
##' 
##' \code{\link{map.identify}} \code{\link{map.sel.poly}}
##' @keywords iplot
##' @rdname spc-identify
##' @export
##' @examples
##' 
##' if (interactive ()){
##' ispc <- sample (nrow (laser), 10)
##' ispc
##' 
##' identified <- spc.identify (plotspc (laser[ispc]))
##' 
##' ## convert to the "real" spectra indices
##' ispc [identified$ispc]
##' identified$wl
##' identified$spc
##' 
##' ## allow the labels to be plotted into the plot margin
##' spc.identify (plotspc (laser[ispc]), ispc = ispc, xpd = NA) 
##' 
##' spc.identify (plotspc (paracetamol, xoffset = 1100,
##'               wl.range = c (600 ~ 1700, 2900 ~ 3150)),
##'               formatter = spc.label.wlonly)
##' 
##' ## looking for minima
##' spc.identify (plot (-paracetamol, wl.reverse = TRUE),
##'               point.fn = spc.point.min, adj = c (1, 0.5))
##' 
##' }
##' 
spc.identify <- function (x, y = NULL, wavelengths = NULL, ispc = NULL,
                          tol.wl = diff (range (x)) / 200, 
                          tol.spc = diff (range (y)) / 50,
                          point.fn = spc.point.max, # function to find the maximum
                          formatter = spc.label.default, # NULL: suppress labels
                          ..., cex = 0.7, adj = c (0, 0.5), srt = 90, # for the label text
                          warn = TRUE){
	
	if (! interactive ())
		stop ("spc.identify works only on interactive graphics devices.")
		
  if (is.list (x)) {
    if (is.null (wavelengths))
      wavelengths <- x$wavelengths
    if (is.null (y))
      y <- x$y
    x <- x$x
  }

  debuglevel <- hy.getOption ("debuglevel")

  if ((length (x) != length (y)) | (length (x) != length (wavelengths)))
    stop ("x, y, and wavelength need to have the same length.")

  if (is.null (ispc))
    ispc <- row (y)
  else
    ispc <- ispc[row(y)]
  
  pts <- data.frame (ispc = rep (NA, 50), wl = NA, spc = NA)
  pos <- 1
  
  while (! is.null (tmp <- locator (n = 1))){
    wl <- approx (x, wavelengths, tmp$x, rule = 2)$y # return wl_min / wl_max for outside pts.

    if (debuglevel == 2L) {
      points (tmp$x, tmp$y, pch = ".", col = "red")
      rect (tmp$x - tol.wl, tmp$y - tol.spc, tmp$x + tol.wl, tmp$y + tol.spc,
            border = "red", col = NA)
    }

    i.window <- wavelengths >= wl - tol.wl & # window to search for the closest spectrum
                wavelengths <= wl + tol.wl &
                y >= tmp$y - tol.spc &
                y <= tmp$y + tol.spc

    if (! any (i.window)){
      if (warn)
        warning ("No spectra in specified window.")
      else
        pos <- pos + 1

      if (debuglevel == 1L) {
        points (tmp$x, tmp$y, pch = ".", col = "red")
        rect (tmp$x - tol.wl, tmp$y - tol.spc, tmp$x + tol.wl, tmp$y + tol.spc,
              border = "red", col = NA)
      }

    } else {
    
      ## find spectrum closest to clicked point.
      ## x and y distances are scaled according to tolerance.
      tmp <- ((wl - wavelengths [i.window]) / tol.wl)^2 +
            ((tmp$y - y [i.window]) / tol.spc)^2 
      tmp <- which (i.window) [which.min (tmp)]

      pts [pos, "ispc"] <- ispc [tmp]   # closest spectrum;
                                        # this will grow the data.frame if necessary
                                        # no time concern with hand-clicked points

      ## search for the max (min) of spectrum pt within tmp$x +- tol.wl
      i.window <- which (ispc == ispc [tmp] &
                         wavelengths >= wl - tol.wl &
                         wavelengths <= wl + tol.wl)
      
      pts [pos, 2 : 3] <- point.fn (wl = wavelengths [i.window],
                                    spc = y [i.window],
                                    wlclick = wl)
      
      ## label the point
      if (! is.null (formatter)){
        lab <- formatter (pts [pos, 1], pts [pos, 2], pts [pos, 3])
        
        text (approx (wavelengths, x, pts [pos, 2], rule = 2),
              pts [pos, 3], labels = lab, cex = cex, adj = adj, srt = srt, ...)
      }

      pos <- pos + 1
    }
    

  }

  pts [seq_len (pos - 1),]
}

##' @rdname spc-identify
##' @param wl the wavelength to label
##' @param spc the intensity to label
##' @param wlclick the clicked wavelength
##' @export
spc.point.max <- function (wl, spc, wlclick){
  i <- which.max (spc)
  c (wl = wl [i], spc = spc [i])
}

##' @rdname spc-identify
##' @export
spc.point.default <- function (wl, spc, wlclick){
  i <- round (approx (wl, seq_along (wl), wlclick, rule = 2)$y)
  c (wl = wl [], spc = spc [i])
}

##' @rdname spc-identify
##' @export
spc.point.min <- function (wl, spc, wlclick){
  i <- which.min (spc)
  c (wl = wl [i], spc = spc [i])
}

##' @rdname spc-identify
##' @export
spc.point.sqr <- function (wl, spc, wlclick, delta = 1L){
  i <- which.max (spc) 

  ## points (wl [i], spc [i])              
  if (i > 1L && i < length (wl)) {
    i <- i + (-delta : delta)
    i <- i %in% seq_along (wl)          # make sure the indices exist
    
    p <- outer (wl [i], 0 : 2, "^")     # Vandermonde matrix
    p <- qr.solve (p, spc [i])

    i <- -p [2] / p [3] / 2

    ## lines (wl, outer (wl, 0 : 2, "^") %*% p, col = "red") 
    c (wl = i, spc = sum (p * c(1, i, i^2)))

  } else {

    c (wl = wl [i], spc = spc [i])

  }
}

##' @param ispc if a selection of spectra was plotted, their indices can be
##'   given in \code{ispc}. In this case \code{ispc [i]} is returned rather
##'   than \code{i}.
##' @param digits how many digits of the wavelength should be displayed?
##' @rdname spc-identify
##' @export
spc.label.default <- function (ispc, wl, spc, digits = 3){
  sprintf(" %i, %s ", ispc, format (wl, digits = digits))
}

##' @rdname spc-identify
##' @export
spc.label.wlonly <- function (ispc, wl, spc, digits = 3){
  sprintf(" %s ", format (wl, digits = digits))
}




