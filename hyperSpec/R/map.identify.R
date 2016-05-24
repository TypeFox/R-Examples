##' @aliases levelplot,hyperSpec,missing-method
##' @include plotmap.R
##' @rdname levelplot
##' @export
##' @seealso  \code{\link[hyperSpec:options]{hyperSpec options}} \code{\link{spc.identify}}
##' \code{\link{map.sel.poly}}
##' @param tol tolerance for \code{map.identify} as fraction of the viewport
##'   (i.e. in "npc" \link[grid]{unit}s)
##' @param warn should a warning be issued if no point is within the specified
##'   tolerance? See also details.
##' @importFrom grid convertX convertY grid.locator grid.circle gpar 
##' @importFrom lattice trellis.focus ltext
map.identify <- function (object, model = spc ~ x * y, voronoi = FALSE, ...,
                          tol = .02, warn = TRUE){
	
	if (! interactive ())
		stop ("map.identify works only on interactive graphics devices.")
	
  chk.hy (object)
  validObject (object)

  dots <- modifyList (list (object = object, model = model, ...),
                      list (subscripts = TRUE))
  
  if (voronoi) { dots <- modifyList (list (col = "black", border = "#00000080"), dots)
    
    ## we need to mix the spectra, otherwise the voronoi plot does not work with 
    ## complete rectangular maps. mix keeps track of the reordering.
    dots$mix <- FALSE
    mix <- sample (nrow (object))
    dots$object <- object [mix]
    lattice <- do.call (plotvoronoi, dots)
    mix <- order (mix)
  } else {
    lattice <- do.call (plotmap, dots)
    mix <- row.seq (object)
  }

  print (lattice)
  trellis.focus ()

  tol = tol^2
  xn <- lattice$panel.args.common$x [mix]
  yn <- lattice$panel.args.common$y [mix]
  x = as.numeric (convertX (unit (xn, "native"), "npc"))
  y = as.numeric (convertY (unit (yn, "native"), "npc"))

  debuglevel <- hy.getOption ("debuglevel")

  res <- numeric (0)
  repeat {
    tmp <- grid.locator (unit = "npc")
    if (debuglevel == 2L)
      grid.circle (tmp[1], tmp[2], sqrt (tol), gp = gpar (col = "red"))

    if (is.null (tmp))
      break
    
    tmp <- as.numeric (tmp)
    d2 <- (x - tmp [1])^2 + (y - tmp [2])^2
    pt <- which.min (d2)
    if (d2 [pt] <= tol) {
      res <- c (res, pt)
      if (debuglevel >= 1L)
        ltext (xn [pt], yn [pt], label = pt)
    } else {
      if (warn) {
        warning ("No point within tolerance (", tol, " = ",
                 convertX (unit (sqrt (tol), "npc"), "native")," x-units or",
                 convertY (unit (sqrt (tol), "npc"), "native")," y-units).")
        if (debuglevel == 1L)
          grid.circle (tmp[1], tmp[2], sqrt (tol), gp = gpar (col = "red"))
      }
    }
  }

  res
}

