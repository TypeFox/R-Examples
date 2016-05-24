vec <- function(x, n = 1, deg = 90) {
  ## x is a unit vector
  ## n is the number of unit vectors with the given degree to generate
  ## deg is the angle, in degrees, between the x and the generated vector(s)

  ## if the user gives more than 180 degrees (in absolute value) make it less than 180

  p <- length(x)  ## number of dimensions
  x <- x / sqrt( sum(x^2) )  ## make sure x is a unit vector

  ### this is for orthogonal vectors
   fu <- function(y, x) {
     y < y / sqrt( sum(y^2) )
     abs( deg / 180 * pi - acos( sum( x * y ) ) )
   }

  runtime <- proc.time()

  if ( abs( deg ) < 0 || abs(deg) > 180 ) {
    deg <- deg %% 180
  }

  if ( deg == 0 ) {
    crit <- 0
    mat <- x  
  } else if (deg == 180) {
    crit <- 180
    mat <-  -x
  } else if (deg > 0 & deg < 180) {
  
    mat <- matrix(nrow = n, ncol = p)
    crit <- numeric(n)
    
    for (i in 1:n) {
      options(warn = -1)
      pa <- nlm(fu, rnorm(p), x = x)
      y <- pa$estimate
      y <- y / sqrt( sum(y^2) )
      ca <- acos( sum(x * y) ) / pi * 180
      while( abs(ca - deg) > 1e-06 ) {
        options(warn = -1)
        pa <- nlm(fu, rnorm(p), x = x)
        y <- pa$estimate
        y <- y / sqrt( sum(y^2) )
        ca <- acos( sum(x * y) ) / pi * 180      
      }
      mat[i,] <- y
      crit[i] <- ca
    }

  }

  runtime <- proc.time() - runtime

  list(runtime = runtime, crit = crit, mat = mat)
}

