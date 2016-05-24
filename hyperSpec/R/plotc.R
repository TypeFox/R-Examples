###-----------------------------------------------------------------------------
###
###  plotc - plot timeseries, concentration, ... 
###
###  C. Beleites



##' Calibration- and Timeseries Plots, Depth-Profiles and the like
##' \code{plotc} plots intensities of a \code{hyperSpec} object over another
##' dimension such as concentration, time, or a spatial coordinate.
##' 
##' If \code{func} is not \code{NULL}, the summary characteristic is calculated
##' first by applying \code{func} with the respective arguments (in
##' \code{func.args}) to each of the spectra. If \code{func} returns more than
##' one value (for each spectrum), the different values end up as different
##' wavelengths.
##' 
##' If the wavelength is not used in the model specification nor in
##' \code{groups}, nor for specifying \code{subsets}, and neither is
##' \code{func} given, then only the first wavelength's intensities are plotted
##' and a warning is issued.
##' 
##' The special column names \code{.rownames} and \code{.wavelength} may be
##' used.
##' 
##' The actual plotting is done by \code{\link[lattice]{xyplot}}.
##' 
##' @param object the \code{hyperSpec} object
##' @param model the lattice model specifying the plot
##' @param func function to compute a summary value from the spectra to be
##'   plotted instead of single intensities
##' @param func.args further arguments to \code{func}
##' @param groups grouping variable, e.g. \code{.wavelength} if intensities of
##'   more than one wavelength should be plotted
##' @param ... further arguments to \code{\link[lattice]{xyplot}}.
##' @author C. Beleites
##' @seealso \code{\link[lattice]{xyplot}}
##' @keywords hplot
##' @export
##' @import graphics
##' @importFrom lattice xyplot
##' @examples
##' 
##' 
##' ## example 1: calibration of fluorescence 
##' plotc (flu) ## gives a warning
##' 
##' plotc (flu, func = mean)
##' plotc (flu, func = range, groups = .wavelength)
##' 
##' plotc (flu[,,450], ylab = expression (I ["450 nm"] / a.u.))
##' 
##' 
##' calibration <- lm (spc ~ c, data = flu[,,450]$.)
##' summary (calibration)
##' plotc (flu [,, 450], type = c("p", "r"))
##' 
##' conc <- list (c = seq (from = 0.04, to = 0.31, by = 0.01))
##' ci <- predict (calibration, newdata = conc, interval = "confidence", level = 0.999) 
##' 
##' panel.ci <-  function (x, y, ...,
##'                        conc, ci.lwr, ci.upr, ci.col = "#606060") {
##'    panel.xyplot (x, y, ...)
##'    panel.lmline (x, y,...)
##'    panel.lines (conc, ci.lwr, col = ci.col)
##'    panel.lines (conc, ci.upr, col = ci.col)
##' }
##' 
##' plotc (flu [,, 450], panel = panel.ci,
##'        conc = conc$c, ci.lwr = ci [, 2], ci.upr = ci [, 3])
##' 
##' ## example 2: time-trace of laser emission modes
##' cols <- c ("black", "blue", "#008000", "red")
##' wl <- i2wl (laser, c(13, 17, 21, 23))
##' 
##' plotspc (laser, axis.args=list (x = list (at = seq (404.5, 405.8, .1))))
##' for (i in seq_along (wl))
##'    abline (v = wl[i], col = cols[i], lwd = 2)
##' 
##' plotc (laser [,, wl], spc ~ t, groups = .wavelength, type = "b",
##'        col = cols)
##' 
plotc <- function (object, model = spc ~ c, groups = NULL,
                     func = NULL, func.args = list (), ...){
  chk.hy (object)
  validObject (object)

  dots <- list (...)

  if (! is.null (func)) 
    object <- do.call (apply, c (list (object, 1, func), func.args))
  
  ## allow to plot against the row number
  object$.row <- row.seq (object)

  groups <- substitute (groups)
  
  ## find out whether the wavelengths are needed individually,
  ## if not, use only the first wavelength and issue a warning
  parsed.formula <- latticeParseFormula (model,
        as.long.df (object [1, , 1, wl.index = TRUE], rownames = TRUE),
        groups = groups, dimension = 2)
  
  use.c <- parsed.formula$right.name
  use.spc <- parsed.formula$left.name

  if (use.spc == "spc" && nwl (object) > 1 && is.null (func) &&
      !any (grepl (".wavelength", c(as.character (model),
                                    as.character (groups),
                                    as.character (dots$subset))))) {
    object <- object [,, 1, wl.index = TRUE]
    warning ("Intensity at first wavelengh only is used.")
  }

  if (is.null (func))
    ylab <- object@label [[use.spc]]
  else {
    ylab <- substitute (func ())
    ylab [[2]] <- object@label [[use.spc]][[1]]
    for (i in seq_along (func.args)){
      if (names (func.args)[[i]] == "")
        ylab [[i + 2]] <- func.args [[i]]
      else
        ylab [[i + 2]] <- bquote (.(x) == .(y),
                                  list (x = names (func.args) [[i]],
                                        y = as.character (func.args [[i]])))
      
      }
    ylab <- as.expression (ylab)
  }
  
  ## set defaults: axis labels, plot style
  dots <- modifyList (list (xlab = object@label [[use.c]],
                            ylab = ylab,
                            pch = 19),
                      dots)

  ## expand the data.frame
  df <- as.long.df (object, rownames = TRUE, wl.factor = TRUE)

  ## if plots should be grouped or conditioned by wavelength,
  ## it is better to have a factor
  if ((! is.null (parsed.formula$condition) &&
       parsed.formula$condition == ".wavelength") ||
      (! is.null (groups) &&
       as.character (groups) == ".wavelength"))
    df$.wavelength <- as.factor (df$.wavelength)

  ## plot
  do.call(xyplot, c (list (x = model, data = df, groups = groups), dots))
}

