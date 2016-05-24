`marksum` <-
function (mippp, R = 10, nx = 30, ny = 30) 
{
    
    dataname <- deparse (substitute(mippp))
    verifyclass(mippp, "ppp")
    if (is.marked(mippp) != TRUE) 
        stop("marksum only implemented for **marked** patterns")
   #generate a grid of dimensions nx ny inside the window of the pattern
    grid <- gridcenters(mippp$window, nx, ny)
          options(warn= -1) # do not warning on trimming points
    grid.ppp <- ppp(x = grid$x, y = grid$y, marks = rep(0, length(grid$x)), 
                         window = mippp$window)[mippp$window]
         options(warn= 0)
    prueba <- superimpose(grid.ppp, mippp)

   ##now count all the points within a distance R of the grid points
    marksum <- markstat(prueba, sum, R = R)
    marksum <- marksum[1:grid.ppp$n] #we want only the sums around the grid points!
    pointsum <- markstat(prueba, length, R = R)
    pointsum <- pointsum[1:grid.ppp$n] #we want only the counts around the grid points!
## correction of excess points in pointsum. We must substract the grid points summed
## in the previous markstat
    minus <- markstat(grid.ppp, length, R = R)
    pointsum <- pointsum - minus
    normalized.marksum <- marksum/pointsum
    normalized.marksum[marksum == 0] <- 0
    result <- list(normalized = normalized.marksum, marksum = marksum,
                     pointsum = pointsum,  minus = minus, grid = grid, grid.ppp=grid.ppp,nx = nx, 
                     ny = ny, R = R, window = mippp$window, dataname=dataname)
   class(result) <- c("ecespa.marksum", class(result))
   return(result)
}

