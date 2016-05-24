`plotdat` <-
function(rho,
                    cyclVar,
                    circle,
                    transf,
                    general,
                    grid)
  ## Author: Rene Locher
  ## Version: 2009-03-13
  ## helper function for plot.rose
  {

    ##----------------------------------------
    ## calculating x & y coordinates for data points

    ## for drawing data
    x.dat <- as.vector(sweep(transf(rho)-transf(grid$ray$lim[1]),
                             MARGIN = 1,
                             sin(2*pi*(cyclVar+general$shift)/circle),
                             "*"))
    y.dat <- as.vector(sweep(transf(rho)-transf(grid$ray$lim[1]),
                             MARGIN = 1,
                             cos(2*pi*(cyclVar+general$shift)/circle),
                             "*"))
    id.dat <- rep(1:ncol(rho), rep(nrow(rho),ncol(rho)))

    return(list(x = x.dat,             ## x-coordinate of point
                y = y.dat,             ## y-coordinate of point
                id = id.dat))           ## category
  } ## plotdat

