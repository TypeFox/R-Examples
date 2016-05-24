# determine points to set to zeros in empty space around some points

boundary <- function(points
                    , density = 0.02
                    , grid = 10
                    , box.offset = 0.1
                    , tightness = "auto"
                    , manual = NULL
                    , plot = TRUE) {

  p <- xy.coords(points)

  if (tightness == "auto") {
    tightness <- MASS::bandwidth.nrd(p$x)
  }

  k <- MASS::kde2d(p$x, p$y, h = tightness, n = grid)

  zeros <- which(k$z < density, arr.ind = TRUE)
  zeroX <- k$x[zeros[,1]]
  zeroY <- k$y[zeros[,2]]

  rX <- diff(range(p$x))*box.offset
  rY <- diff(range(p$y))*box.offset
  pXmin <- min(k$x)-rX
  pXmax <- max(k$x)+rX
  pYmin <- min(k$y)-rY
  pYmax <- max(k$y)+rY

  borderX <- c(  k$x
               , k$x
               , rep(pXmin, times = length(k$y))
               , rep(pXmax, times = length(k$y))
               , pXmin, pXmin, pXmax, pXmax
               )
  borderY <- c(  rep(pYmin, length(k$x))
               , rep(pYmax, length(k$x))
               , k$y
               , k$y
               , pYmin, pYmax, pYmin, pYmax
               )

  if (plot) {

    plot(borderX, borderY, col = "blue", pch = 19, xlab = "", ylab = "")
    title(xlab = paste( "density =", density
                      , "grid =", grid
                      , "\nbox.offset = ", box.offset
                      , "tighness = ", round(tightness, 1)
                      ), col = "grey", cex = 0.5)
    contour(k, add = TRUE)
    contour(k, level = density, col = "red", add = TRUE)
    points(zeroX, zeroY, col = "red", pch =  19)
    points(manual, col = "green", pch = 19)
    points(points, pch = 20)

  } else {

    return(cbind( x = c(zeroX, borderX, manual[,1])
                , y = c(zeroY, borderY, manual[,2])
                ))
  }
}
