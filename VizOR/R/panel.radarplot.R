## (Adapted from panel.xyplot)
##' Panel function for multidimensional 'Table 1 at-a-glance'
##' See \code{\link{radarplot}}
##' 
##' @param x
##' @param y
##' @param groups
##' @param pch
##' @param col
##' @param col.line
##' @param col.symbol
##' @param font
##' @param fontfamily
##' @param fontface
##' @param lty
##' @param cex
##' @param fill
##' @param lwd
##' @param \dots
##' @author David C. Norris
##' @seealso \code{\link{radarplot}}
##' @keywords internal hplot
##' @export panel.radarplot
panel.radarplot <-
  function (x, y, groups = NULL, pch = if (is.null(groups)) plot.symbol$pch else superpose.symbol$pch, 
            col, col.line = if (is.null(groups)) plot.line$col else superpose.line$col, 
            col.symbol = if (is.null(groups)) plot.symbol$col else superpose.symbol$col, 
            font = if (is.null(groups)) plot.symbol$font else superpose.symbol$font, 
            fontfamily = (if (is.null(groups)) plot.symbol else superpose.symbol)$fontfamily, 
            fontface = (if (is.null(groups)) plot.symbol else superpose.symbol)$fontface, 
            lty = if (is.null(groups)) plot.line$lty else superpose.line$lty, 
            cex = if (is.null(groups)) plot.symbol$cex else superpose.symbol$cex, 
            fill = if (is.null(groups)) plot.symbol$fill else superpose.symbol$fill, 
            lwd = if (is.null(groups)) plot.line$lwd else superpose.line$lwd,
            w = 0.12, h = 0.06, # TODO: Eliminate these 'tweaks' by calculating text size directly
            strength=NULL,
            ...)
{
  if (all(is.na(x) | is.na(y))) 
    return()
  plot.symbol <- trellis.par.get("plot.symbol")
  plot.line <- trellis.par.get("plot.line")
  superpose.symbol <- trellis.par.get("superpose.symbol")
  superpose.line <- trellis.par.get("superpose.line")
  if (!missing(col)) {
    if (missing(col.line)) 
      col.line <- col
    if (missing(col.symbol)) 
      col.symbol <- col
  }
  n <- length(list(...)$radii)
  theta <- (pi/2)*(1 - 2/n) + (0:(n-1))*(2*pi)/n
  if (!is.null(groups)) {
    panel.segments(x0 = 0, y0 = 0, x1 = cos(theta), y1 = sin(theta),
                   col = 'gray', lwd = lwd, ...)
    ## If 'strength' argument supplied in '...' then
    if(!is.null(strength)){
      ## Let the spokes be drawn as sectors with angle proportional to strength
      max.angle <- 0.1 * (2*pi)/length(x) # let strongest association spread to 10% of sector
      phi <- max.angle*strength/max(strength)
      fan <- (phi/2) %o% seq(-1,1,0.1) # matrix of segment arc displacements
      xs <- as.vector(t(cbind(0, cos(theta+fan))))
      ys <- as.vector(t(cbind(0, sin(theta+fan))))
      grid.polygon(x = unit(xs, "native"),
                   y = unit(ys, "native"),
                   gp = gpar(col="gray", fill="gray"))
      ## Then plot strength-of-association on radii
      if(is.null(names(strength)))
        names(strength) <- format(round(strength, digits=2), digits=2)
      grid.text(label = names(strength),
                x = unit((1 + w)*cos(theta), "native"),
                y = unit((1 + h)*sin(theta), "native"))
    }
    grid.text(label=list(...)$radii,
              x = unit(0.5*cos(theta),"native"),
              y = unit(0.5*sin(theta),"native"),
              rot = ((180/pi)*theta + 90) %% 180 - 90,
              just = "center")
    panel.superpose(x, y, groups = groups, pch = pch, 
                    col.line = col.line, col.symbol = col.symbol, font = font, 
                    fontfamily = fontfamily, fontface = fontface, lty = lty, 
                    cex = cex, fill = fill, lwd = lwd,
                    panel.groups = panel.radarplot,
                    grid = FALSE, ...)
  }
  else {
    r <- as.numeric(y)
    x <- r * cos(theta)
    y <- r * sin(theta)
    panel.polygon(x = x, y = y, lty = lty, border = col.line, lwd = lwd, ...)
  }
}
