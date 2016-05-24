shadeDens = function(x0, x1, dens, col = "lightgrey", n.points  = 200,
                     lty = 1,...){
    xPts = seq(x0, x1, length = n.points)
    yPts = dens(xPts,...)
    xPts = c(x0, xPts, x1)
    yPts = c(0, yPts, 0)
    polygon(xPts, yPts, col = col, lty = lty)
}
