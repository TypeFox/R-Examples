#' @title Pizza color wheel
#' 
#' @description
#' This function displays a color wheel with specified colors
#' 
#' @details
#' This function is based on the \code{\link{pie}} function
#' 
#' @param colors a vector with R color names of colors in hexadecimal notation
#' @param bg background color of the plot. Default \code{"gray95"}
#' @param border color of the border separating the pizza slices
#' @param init.angle integer value indicating the start angle (in degrees) for
#' the slices
#' @param cex numeric value indicating the character expansion of the labels
#' @param lty argument passed to \code{\link{polygon}} which draws each slice
#' @param labcol color for the labels (i.e. names of the colors)
#' @param \dots graphical parameters (\code{\link{par}}) can be given as
#' argument to \code{pizza}
#' @author Gaston Sanchez
#' @seealso \code{\link{wheel}}
#' @export
#' @examples
#'
#' # pizza color wheel for rainbow colors
#' pizza(rainbow(7))
#' 
#' # pizza color wheel for tomato (18 colors)
#' pizza(setColors("tomato", 18), bg = "gray20", cex = 0.7)
#'
pizza <-
function(colors, bg = "gray95", border = NA, 
	init.angle = 105, cex = 0.8, lty = 1, labcol = NULL, ...)
{
  n <- length(colors)
  x <- rep(1, n)
  x <- c(0, cumsum(x)/sum(x))
  dx <- diff(x)
  nx <- length(dx)
  # set colors
  labels = colors
  # 
  if (is.null(labcol))
  {
    if (mean(col2rgb(bg)) > 127)
      labcol = rep("black", n)
    if (mean(col2rgb(bg)) <= 127)
      labcol = rep("white", n)
  }
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
  border <- rep(border, length.out = nx)
  if (is.null(border[1]))
    border <- rep(bg, length.out = nx)
  lty <- rep(lty, length.out = nx)
  angle <- rep(45, length.out = nx)
  radius = seq(1, 0, by=-1/n)[1:n]
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
            border = border[i], col = colors[i], lty = lty[i])
    P <- t2xy(mean(x[i + 0:1]), rad=radius[1])
    lab <- labels[i]
    if (!is.na(lab) && nzchar(lab)) {
      adjs = 0.5
      if (P$x > 1e-08) adjs <- 0
      if (P$x < -1e-08) adjs <- 1
      lines(c(1, 1.05) * P$x, c(1, 1.05) * P$y, col=labcol[i])
      text(1.1 * P$x, 1.1 * P$y, labels[i], xpd = TRUE, 
           adj = adjs, cex=cex, col=labcol[i], ...)
    }
  }
  invisible(NULL)
}
