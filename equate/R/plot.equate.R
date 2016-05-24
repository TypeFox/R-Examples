#' Plotting Equating Results
#' 
#' Functions for plotting equating functions from one or more objects
#' of class \dQuote{\code{equate}} or \dQuote{\code{equate.list}}.
#' 
#' Equating functions (\code{out = "eqs"}) are plotted as lines based on
#' the concordance table for each equating object that is supplied. Standard
#' errors (\code{out = "se"}) default to bootstrap standard errors, if
#' available, otherwise, analyitical standard errors are plotted. Bias
#' (\code{out = "bias"}) and RMSE (\code{out = "rmse"}) are also taken
#' from bootstrapping output.
#' 
#' @param ... one or more equating objects, each containing results for
#' equating the same two test forms.
#' @param elist list of equatings to be plotted.
#' @param add logical, with default \code{FALSE}, specifying whether to
#' create a new plot or add to the current one.
#' @param out character vector specifying the output to be plotted, either
#' equating functions (\code{"eqs"}), standard errors (\code{"se"}),
#' bias (\code{"bias"}), or RMSE (\code{"rmse"}).
#' @param xpoints,ypoints optional vectors of the same length containing
#' raw scores on forms X and Y, assuming a single group or equivalent groups
#' design.
#' @param addident logical, with default \code{TRUE}, for plotting the
#' identity function. The result depends on \code{out}.
#' @param identy vector of y coordinates for plotting the identity line.
#' Defaults to the X scale when \code{out = "eqs"}, otherwise, a horizontal
#' line with intercept 0.
#' @param identcol color used for plotting the identity line.
#' @param rescale intercept and slope, with default 0 and 1, used to rescale
#' all lines before plotting.
#' @param xlab,ylab,col,pch,lty,lwd graphical parameters passed to \code{par},
#' with \code{col}, \code{pch}, \code{lty}, and \code{lwd} recycled as necessary.
#' @param subset vector for subsetting the output when multiple equating
#' functions are included in \code{x}.
#' @param morepars list of additional graphical parameters, excluding
#' \code{xlab}, \code{ylab}, \code{col}, \code{pch}, \code{lty}, and \code{lwd}.
#' @param addlegend logical, with default \code{TRUE}, indicating whether or
#' not a legend should be added.
#' @param legendtext character vector of text to be passed to the \code{legend}
#' argument of the \code{legend} function, defaulting to a combination of the
#' equating types and methods specified in each equating object.
#' @param legendplace placement of the legend.
#' @param x \dQuote{\code{\link{equate.list}}} object, containing output from
#' multiple equatings.
#' @examples
#' 
#' # See ?equate for additional examples
#' 
#' rx <- as.freqtab(ACTmath[, 1:2])
#' ry <- as.freqtab(ACTmath[, c(1, 3)])
#' set.seed(2007)
#' 
#' req1 <- equate(rx, ry, type = "i", boot = TRUE, reps = 5)
#' req2 <- equate(rx, ry, type = "m", boot = TRUE, reps = 5)
#' req3 <- equate(rx, ry, type = "l", boot = TRUE, reps = 5)
#' req4 <- equate(rx, ry, type = "e", boot = TRUE, reps = 5,
#'   smooth = "loglin", degree = 3)
#' req5 <- composite(list(req1, req2), wc = .5, symmetric = TRUE)
#' 
#' plot(req1, req2, req3, req4, req5[[1]], addident = FALSE)
#' plot(req5)
#' 
#' @export
plot.equate <- function(..., elist = NULL, add = FALSE,
  out = "eqs", xpoints, ypoints, addident = TRUE,
  identy, identcol = 1, rescale = c(0, 1),
  xlab = "Total Score", ylab, col = rainbow(length(x)),
  pch, lty = 1, lwd = 1, subset, morepars = NULL, addlegend = TRUE,
  legendtext, legendplace = "bottomright") {
  
  x <- c(list(...), elist)
  if(missing(subset)) subset <- 1:length(x)
  x <- x[subset]
  nx <- length(x)
  xscale <- scales(x[[1]]$x)
  
  out <- match.arg(tolower(out),
    c("se", "bias", "eqs", "rmse"))
  if(out == "se") {
    y <- lapply(x, function(z) {
      if(is.null(z$bootstraps)) z$se
      else z$bootstraps$se
    })
  }
  else if(out == "bias")
    y <- lapply(x, function(z) z$bootstraps$bias)
  else if(out == "rmse")
    y <- lapply(x, function(z) z$bootstraps$rmse)
  else if(out == "eqs")
    y <- lapply(x, function(z) z$concordance[, 2])
  else
    stop("'out' must be one of 'eqs', 'se', 'bias' ",
      "or rmse")
  if(any(unlist(lapply(y, is.null))))
    stop("one or more equatings does not contain ", out)
  y <- lapply(y, function(z) z*rescale[2] + rescale[1])
  
  if(missing(ylab))
    ylab <- switch(out, eqs = "Equated Score",
      se = "Standard Error", bias = "Bias",
      rmse = "RMSE")
  if(!is.null(morepars)) {
    nopars <- c("xlab", "ylab", "col", "lty", "pch", "lwd")
    noparsl <- nopars %in% names(morepars)
    if(any(noparsl)) {
      warning("the following graphical parameter(s)",
        " must be specified outside of 'morepars': ",
        paste(nopars[noparsl], collapse = ", "))
      morepars <- morepars[!names(morepars) %in%
          nopars]
    }
  }
  
  if(!add) {
    do.call(plot, c(list(x = range(xscale),
      y = range(y), xlab = xlab, ylab = ylab,
      type = "n"), morepars))
    if(!missing(xpoints) && is.freqtab(xpoints))
      do.call(points.freqtab, c(list(x = xpoints,
        xcol = "lightgray"), morepars))
    else if(!missing(xpoints) & !missing(ypoints))
      do.call(points, c(list(x = xpoints,
        y = ypoints, col = "lightgray"), morepars))
  }
  
  if(addident) {
    if(missing(identy))
      identy <- switch(out, eqs = xscale,
        rep(0, length(xscale)))
    lines(xscale, identy*rescale[2] + rescale[1],
      col = identcol)
  }
  
  col <- rep(col, length = nx)
  lty <- rep(lty, length = nx)
  lwd <- rep(lwd, length = nx)
  for(i in 1:nx)
    lines(xscale, y[[i]], col = col[i],
      lty = lty[i], lwd = lwd[i])
  if(!missing(pch)) {
    pch <- rep(pch, length = nx)
    for(i in 1:nx)
      points(xscale, y[[i]], col = col[i],
        pch = pch[i])
  }
  
  if(addlegend) {
    if(missing(legendtext)) {
      legendtext <- c("identity", "mean", "linear",
        "general", "circle", "equip", "composite")
      legendtext <- lapply(x, function(z)
        legendtext[charmatch(substr(z$type, 1, 2),
          legendtext)])
      if(x[[1]]$design == "nonequivalent groups") {
        methods <- c("nW", "chain", "b/H", "tucker",
          "levine", "fE")
        methods <- lapply(x, function(z)
          methods[charmatch(substr(z$method, 1, 1),
            methods)])
        legendtext <- paste(legendtext,
          methods, sep = ": ")
        legendtext[grep("ident", legendtext)] <-
          "identity"
        legendtext[grep("comp", legendtext)] <-
          "composite"
      }
      legendtext <- gsub("\\b(\\w)", "\\U\\1",
        legendtext, perl = TRUE)
    }
    if(addident) {
      legendtext <- c("Identity", legendtext)
      lty = c(1, lty)
      col = c(identcol, col)
    }
    legend(legendplace, legend = legendtext,
      lty = lty, col = col, bty = "n")
  }
}

# @describeIn plot.equate Method for plotting \dQuote{\code{equate.list}}
# objects directly.
#' @rdname plot.equate
#' @export
plot.equate.list <- function(x, ...) {
  
  plot.equate(elist = x, ...)
}