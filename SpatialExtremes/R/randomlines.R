vdc <- function(n, rand.rot = FALSE){
  lines <- .C("vandercorput", as.integer(n), coord = double(3 * n),
              PACKAGE = "SpatialExtremes")$coord

  if (rand.rot){
    ##Gerenerate a random rotation
    u <- runif(3, -.5, .5)
    u <- u / sqrt(sum(u^2))
    angle <- runif(1, 0, 2 * pi)
        
    lines <- .C("rotation", lines = as.double(lines), as.integer(n),
                as.double(u[1]), as.double(u[2]), as.double(u[3]),
                as.double(angle), PACKAGE = "SpatialExtremes")$lines
  }

  lines <- matrix(lines, ncol = 3)
  return(lines)
}
