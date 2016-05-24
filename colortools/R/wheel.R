#' @title Color Wheel
#' 
#' @description
#' This function generates a color wheel for a given color
#' 
#' @details
#' This function is based on the \code{\link{pie}} function
#' 
#' @param color an R color name or a color in hexadecimal notation
#' @param num integer value indicating how many colors to be generated for the
#' color wheel
#' @param bg background color of the plor
#' @param border color of the border separating the slices
#' @param init.angle integer value indicating the start angle (in degrees) for
#' the slices
#' @param cex numeric value indicating the character expansion of the labels
#' @param lty argument passed to \code{\link{polygon}} which draws the slices
#' @param main an overall title for the plot
#' @param verbose logical value indicating whether to return the color names
#' @param \dots graphical parameters (\code{\link{par}}) can be given as
#' argument to \code{wheel}
#' @return A character vector with the color names of the generated wheel in
#' hexadecimal notation
#' @author Gaston Sanchez
#' @seealso \code{\link{pizza}}
#' @export
#' @examples
#' # wheel color with 18 slices for 'tomato'
#' wheel("tomato", num = 18, bg = "gray20", cex = 0.7)
#'
wheel <-
function(color, num=12, bg="gray95", border=NULL, 
	init.angle=105, cex=1, lty=NULL, main=NULL, verbose=TRUE, ...)
{
  if (!is.numeric(num) || any(is.na(num) | num < 0)) 
    stop("\n'num' must be positive")
  x <- rep(1, num)
  x <- c(0, cumsum(x)/sum(x))
  dx <- diff(x)
  nx <- length(dx)
  # set colors
  col = setColors(color, num)
  labels = col
  # labels color
  labcol = ifelse( mean(col2rgb(bg)) > 127, "black", "white")
  # prepare plot window
  par(bg = bg)
  plot.new()
  pin <- par("pin")
  xlim <- ylim <- c(-1, 1)
  if (pin[1L] > pin[2L]) 
    xlim <- (pin[1L]/pin[2L]) * xlim
  else ylim <- (pin[2L]/pin[1L]) * ylim
  dev.hold()
  on.exit(dev.flush())
  plot.window(xlim, ylim, "", asp = 1)
  # get ready to plot
  if (is.null(border[1])) {
    border <- rep(bg, length.out = nx)    
  } else {
    border <- rep(border, length.out = nx)    
  }
  if (!is.null(lty))
    lty <- rep(NULL, length.out = nx)
  angle <- rep(45, length.out = nx)
  radius = seq(1, 0, by=-1/num)[1:num]
  twopi <- -2 * pi
  t2xy <- function(t, rad) {
    t2p <- twopi * t + init.angle * pi/180
    list(x = rad * cos(t2p), y = rad * sin(t2p))
  }
  # plot colored segments
  for (i in 1L:nx)
  {
    n <- max(2, floor(200 * dx[i]))
    P <- t2xy(seq.int(x[i], x[i + 1], length.out = n), rad=radius[1])
    polygon(c(P$x, 0), c(P$y, 0), angle = angle[i], 
            border = border[i], col = col[i], lty = lty[i])
    P <- t2xy(mean(x[i + 0:1]), rad=radius[1])
    lab <- labels[i]
    if (!is.na(lab) && nzchar(lab)) {
      adjs = 0.5
      if (P$x > 1e-08) adjs <- 0
      if (P$x < -1e-08) adjs <- 1
      lines(c(1, 1.05) * P$x, c(1, 1.05) * P$y)
      text(1.1 * P$x, 1.1 * P$y, labels[i], xpd = TRUE, 
           adj = adjs, cex=cex, col=labcol, ...)
    }
  }
  # add title
  title(main = main, ...)
  # return color names
  if (verbose)
    col
}
