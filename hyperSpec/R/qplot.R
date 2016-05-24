##' Spectra plotting with ggplot2
##'
##' These functions are still experimental and may change substantially in future.
##' @title Spectra plotting with ggplot2
##' @param x hyperSpec object
##' @param wl.range wavelength ranges to plot
##' @param ... handed to \code{\link[ggplot2]{geom_line}}
##' @param mapping see  \code{\link[ggplot2]{geom_line}}
##' @param spc.nmax maximum number of spectra to plot
##' @param map.lineonly if \code{TRUE}, \code{mapping} will be handed to
##' \code{\link[ggplot2]{geom_line}} instead of \code{\link[ggplot2]{ggplot}}.
##' @return a \code{\link[ggplot2]{ggplot}} object
##' @author Claudia Beleites
##' @export 
##' @seealso \code{\link{plotspc}}
##'
##' \code{\link[ggplot2]{ggplot}}\code{\link[ggplot2]{geom_line}}
##' @examples
##' 
##'   qplotspc (chondro)
##'   qplotspc (paracetamol, c (2800 ~ max, min ~ 1800)) + scale_x_reverse (breaks = seq (0, 3200, 400))
##' 
##'   qplotspc (aggregate (chondro, chondro$clusters, mean),
##'             mapping = aes (x = .wavelength, y = spc, colour = clusters)) +
##'     facet_grid (clusters ~ .)
##' 
##'   qplotspc (aggregate (chondro, chondro$clusters, mean_pm_sd),
##'             mapping = aes (x = .wavelength, y = spc, colour = clusters, group = .rownames)) +
##'     facet_grid (clusters ~ .)
qplotspc <- function (x,
                      wl.range, ...,
                      mapping = aes_string (x = ".wavelength", y = "spc", group = ".rownames"),
                      spc.nmax = 10,
                      map.lineonly = FALSE){
  chk.hy (x)
  validObject (x)
  
  ## cut away everything that isn't asked for _before_ transforming to data.frame
  if (nrow (x) > spc.nmax) {
    warning ("Number of spectra exceeds spc.nmax. Only the first ", spc.nmax, " are plotted.")
    x <- x [seq_len (spc.nmax)]
  }
  
  if (!missing (wl.range))
    x <- x [,, wl.range]

  df <- as.long.df (x, rownames = TRUE, na.rm = FALSE) # with na.rm trouble with wl.range

  # different spectral ranges go into facets
  if (!missing (wl.range)){
    ranges <- integer (nwl (x))
    for (r in seq_along (wl.range)) 
      ranges [wl2i (x, wl.range [r])] <- r
    if (any (ranges == 0))
      stop ("internal error in qplotspc: 0 range. Please contact the package maintainer.")
    ranges <- as.factor (ranges)
    df$.wl.range <- rep (ranges, each = nrow (x))
  }

  df <- df [! is.na (df$spc),, drop = FALSE]
  if (map.lineonly)
      p <- ggplot (df) + geom_line (mapping = mapping, ...)
  else
      p <- ggplot (df, mapping = mapping) + geom_line (...)

  p <- p + xlab (labels (x, ".wavelength")) + ylab (labels (x, "spc")) 

  if (!missing (wl.range))
      p <- p + facet_grid (. ~ .wl.range, labeller = function (...) "",
                           scales = "free", space = "free") +
               theme (strip.text.x = element_text (size = 0))

  p
}


##' Spectra plotting with ggplot2
##'
##' These functions are still experimental and may change substantially in future.
##'
##' Note that \code{qplotmap} will currently produce the wrong scales if x or y are discrete. 
##' @title Spectra plotting with ggplot2
##' @param object  hyperSpec object
##' @param mapping see  \code{\link[ggplot2]{geom_tile}}
##' @param ... handed to \code{\link[ggplot2]{geom_tile}}
##' @param func function to summarize the wavelengths
##' @param func.args arguments to \code{func}
##' @param map.tileonly if \code{TRUE}, \code{mapping} will be handed to
##' \code{\link[ggplot2]{geom_tile}} instead of \code{\link[ggplot2]{ggplot}}.
##' @return a \code{\link[ggplot2]{ggplot}} object
##' @export 
##' @author Claudia Beleites
##' @seealso \code{\link{plotmap}}
##'
##' \code{\link[ggplot2]{ggplot}}\code{\link[ggplot2]{geom_tile}}
##' @examples
##' qplotmap (chondro)
##' qplotmap (chondro) + scale_fill_gradientn (colours = alois.palette ())
##'
##' ## works also with discrete x or y axis:
##' qplotmap (chondro, mapping = aes (x = x, y = as.factor (y), fill = spc)) 
qplotmap <- function (object, mapping = aes_string (x = "x", y = "y", fill = "spc"), ...,
                      func = mean, func.args = list (),
                      map.tileonly = FALSE){
  chk.hy (object)
  validObject (object)

  if (nwl (object) > 1 & ! is.null (func))
    object <- do.call (apply, c (list (object, 1, func), func.args))

  if (map.tileonly)
      p <- ggplot (as.long.df (object)) + geom_tile (mapping = mapping) 
  else
      p <- ggplot (as.long.df (object), mapping = mapping) + geom_tile () 
  
  p <- p + coord_equal ()

  ## set expand to c(0, 0) to suppress the gray backgroud
  if (is.factor (with (p$data, eval (p$mapping$x))))
        p <- p + scale_x_discrete (expand = c(0, 0)) 
  else
        p <- p + scale_x_continuous (expand = c(0, 0)) 

  if (is.factor (with (p$data, eval (p$mapping$y))))
        p <- p + scale_y_discrete (expand = c(0, 0)) 
  else
      p <- p + scale_y_continuous (expand = c(0, 0)) 

  ## generate axis/scale labels
  ## TODO: own function
  x <- as.character (mapping$x)
  xlabel <- labels (object)[[tail (x, 1)]]
  if (is.null (xlabel)) xlabel <- x

  y <- as.character (mapping$y)
  ylabel <- labels (object)[[tail (y, 1)]]
  if (is.null (ylabel)) ylabel <- y

  f <- as.character (mapping$fill)
  flabel <- labels (object)[[tail (f, 1)]]
  if (is.null (flabel)) flabel <- f

  p + labs (x = xlabel, y = ylabel, fill = flabel) 
}


##' Spectra plotting with ggplot2
##'
##' These functions are still experimental and may change substantially in future.
##' @title Spectra plotting with ggplot2
##' @param object hyperSpec object
##' @param mapping see  \code{\link[ggplot2]{geom_point}}
##' @param ... handed to \code{\link[ggplot2]{geom_point}}
##' @export 
##' @param func function to summarize the wavelengths, if \code{NULL}, only the first wavelength is used
##' @param func.args arguments to \code{func}
##' @param map.pointonly if \code{TRUE}, \code{mapping} will be handed to
##' \code{\link[ggplot2]{geom_point}} instead of \code{\link[ggplot2]{ggplot}}.
##' @return a \code{\link[ggplot2]{ggplot}} object
##' @author Claudia Beleites
##' @seealso \code{\link{plotc}}
##'
##' \code{\link[ggplot2]{ggplot}}\code{\link[ggplot2]{geom_point}}
##' @examples
##' qplotc (flu)
##' qplotc (flu) + geom_smooth (method = "lm")
qplotc <- function (object, mapping = aes_string(x = "c", y = "spc"), ...,
                    func = NULL, func.args = list (),
                    map.pointonly = FALSE){
  chk.hy (object)
  validObject (object)

  dots <- list (...)

  if (! is.null (func)) 
    object <- do.call (apply, c (list (object, 1, func), func.args))
  
  ## allow to plot against the row number
  object$.row <- seq (object, index = TRUE)

  ## find out whether the wavelengths are needed individually,
  ## if not, use only the first wavelength and issue a warning

  if (any (grepl ("spc", as.character (mapping))) && # use intensities
      nwl (object) > 1 &&                            # has > 1 wavelength
      is.null (func) &&                              # no stats function
      ! any (grepl ("[.]wavelength", as.character (mapping)))) {
    object <- object [,, 1, wl.index = TRUE]
    warning ("Intensity at first wavelengh only is used.")
  }

  ## produce fancy y label
  ylab <- labels (object, as.character (mapping$y))
  if (! is.null (func)) 
    ylab <- make.fn.expr (substitute (func), c (ylab, func.args))
  ylab <- as.expression (ylab)
  
  ## expand the data.frame
  df <- as.long.df (object, rownames = TRUE, wl.factor = TRUE)

  ## if plots should be grouped, faceted, etc. by wavelength, it is better to have a factor
  if (any (grepl ("[.]wavelength", mapping [! names (mapping) %in% c("x", "y")])))
    df$.wavelength <- as.factor (df$.wavelength)

  if (map.pointonly)
      p <- ggplot (df) + geom_point (mapping = mapping)
  else
      p <- ggplot (df, mapping = mapping) + geom_point ()
  
  p + ylab (ylab) +
      xlab (labels (object, as.character (mapping$x)))
}

make.fn.expr <- function (fn, l = list ()){

  if (length (fn) > 1L)
    fn <- "f"

  l <- lapply (l, function (x) if (is.logical (x)) as.character (x) else x)

  if (is.null (names (l)))
    names (l) <- rep ("", length (l))
  
  tmp <- mapply (function (x, y) if (nzchar (x) > 0L) bquote (.(x) == .(y)) else y,
                 names (l), l)
  
  e <- expression (f (x))
  e [[1]][[1]] <- fn
  if (length (tmp) > 0L)
    e [[1]][seq_along (tmp) + 1] <- tmp
  else
    e [[1]][2] <- NULL

  e
}
