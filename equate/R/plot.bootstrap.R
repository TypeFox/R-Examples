#' Plotting Bootstrap Equating Results
#' 
#' This function plots bootstrap equating results for objects of class
#' \dQuote{\code{\link{bootstrap}}}.
#' 
#' Lines are plotted for the chosen output type, whether mean equated scores
#' across replications (\code{out = "mean"}), standard errors (\code{out =
#' "se"}), bias (\code{out = "bias"}) or RMSE (\code{out = "rmse"}). The result
#' is similar to that of \code{\link{plot.equate}}.
#' 
#' @param x output from the \code{\link{bootstrap}} function.
#' @param add logical, with default \code{FALSE}, specifying whether to create
#' a new plot or add to the current one.
#' @param out character vector specifying the output to be plotted, either the
#' mean equated scores (\code{"mean"}), standard errors (\code{"se"}), bias
#' (\code{"bias"}), or RMSE (\code{"rmse"}).
#' @param xpoints,ypoints optional vectors of the same length containing raw
#' scores on forms X and Y, assuming a single group or equivalent groups
#' design.
#' @param addident logical, with default \code{TRUE}, for plotting the identity
#' function. The result depends on \code{out}.
#' @param identy vector of y coordinates for plotting the identity line.
#' Defaults to the identity function when \code{out = "eqs"}, otherwise, a
#' horizontal line with intercept 0.
#' @param identcol color used for plotting the identity line.
#' @param rescale intercept and slope, with default 0 and 1, used to rescale
#' all lines before plotting.
#' @param xlab,ylab,col,pch,lty graphical parameters passed to \code{par}, with
#' the lengths of col and lty recycled as necessary.
#' @param subset vector for subsetting the output when multiple equating
#' functions are included in \code{x}.
#' @param morepars list of additional graphical parameters, excluding
#' \code{xlab}, \code{ylab}, \code{col}, \code{pch}, \code{lty}.
#' @param addlegend logical, with default \code{TRUE}, indicating whether or
#' not a legend should be added.
#' @param legendtext character vector of text to be passed to the \code{legend}
#' argument of the \code{legend} function, defaulting to a combination of the
#' equating types and methods specified in each equating object.
#' @param legendplace placement of the legend.
#' @param \dots further arguments passed to or from other methods, excluding
#' graphical parameters.
#' @author Anthony Albano \email{tony.d.albano@@gmail.com}
#' @seealso \code{\link{bootstrap}}, \code{\link{plot.equate}}
#' @keywords misc
#' @examples
#' 
#' set.seed(122713)
#' neat.x <- freqtab(KBneat$x, scales = list(0:36, 0:12))
#' neat.y <- freqtab(KBneat$y, scales = list(0:36, 0:12))
#' eqargs <- list(m.t = list(type = "mean", method = "t"),
#'   l.t = list(type = "lin", method = "t"),
#' 	c.t = list(type = "circ", method = "t"))
#' bootout <- bootstrap(x = neat.x, y = neat.y, args = eqargs,
#' 	reps = 20)
#' plot(bootout, out = "se", legendplace = "top")
#' 
#' @export
plot.bootstrap <- function(x, add = FALSE, out = "mean",
  xpoints, ypoints, addident = TRUE, identy,
  identcol = 1, rescale = c(0, 1), xlab = "Total Score",
  ylab, col = rainbow(length(x$args)), pch, lty = 1,
  subset, morepars = NULL, addlegend = TRUE,
  legendtext, legendplace = "bottomright", ...) {
  
  if(missing(subset)) subset <- 1:length(x$args)
  x$args <- x$args[subset]
  nx <- length(subset)
  xscale <- scales(x$x)
  
  out <- match.arg(tolower(out),
    c("se", "bias", "mean", "rmse"))
  if(!out %in% names(x))
    stop(paste("'x' does not contain", out))
  y <- cbind(cbind(x[[out]])[, subset])
  y <- apply(y, 2, function(z) z*rescale[2] + rescale[1])
  
  if(missing(ylab))
    ylab <- switch(out, mean = "Mean Equated Score",
      se = "Standard Error", bias = "Bias",
      rmse = "RMSE")
  if(!is.null(morepars)) {
    nopars <- c("xlab", "ylab", "col", "lty", "pch")
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
    if(!missing(xpoints) & !missing(ypoints))
      do.call(points, c(list(x = xpoints,
        y = ypoints, col = "lightgray"), morepars))
  }
  
  if(addident) {
    if(missing(identy))
      identy <- switch(out, mean = xscale,
        rep(0, length(xscale)))
    lines(xscale, identy*rescale[2] + rescale[1],
      col = identcol)
  }
  
  lty <- rep(lty, length = nx)
  col <- rep(col, length = nx)
  for(i in 1:nx)
    lines(xscale, y[, i], col = col[i],
      lty = lty[i])
  if(!missing(pch)) {
    pch <- rep(pch, length = nx)
    for(i in 1:nx)
      points(xscale, y[, i], col = col[i],
        pch = pch[i])
  }
  
  if(addlegend) {
    if(missing(legendtext)) {
      legendtext <- c("identity", "mean", "linear",
        "general", "circle", "equip", "composite")
      legendtext <- lapply(x$args, function(z)
        legendtext[charmatch(substr(z$type, 1, 2),
          legendtext)])
      if(margins(x$x) == 2 & !is.null(x$y)) {
        methods <- c("nW", "chain", "b/H", "tucker",
          "levine", "fE")
        methods <- lapply(x$args, function(z)
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
