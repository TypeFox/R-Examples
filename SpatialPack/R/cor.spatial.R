cor.spatial <-
function(x, y, coords)
{
  ## validating arguments
  if (length(x) != length(y))
    stop("'x' and 'y' must have the same length")
  if (!is.numeric(x)) stop("'x' must be a numeric vector")
  if (!is.numeric(y)) stop("'y' must be a numeric vector")
  ## in order to remove all NAs
  OK <- complete.cases(x, y)
  x <- x[OK]
  y <- y[OK]
  n <- length(OK)
  rk.x <- rank(x, ties.method = "first")
  rk.y <- rank(y, ties.method = "first")

  ## is assumed that the coordinates are in the appropiate order
  coords <- as.matrix(coords)
  p <- ncol(coords)
  if (p < 2) stop("'coords' must be a matrix with two columns")
  if (p > 2) warning("only the first two columns of 'coords' are considered")
  p <- 2 # only implemented for this case!
  xpos <- coords[,1]
  ypos <- coords[,2]

  ## initial computations
  dims <- c(n, p, 0)
  bars <- apply(coords, 2, mean)
  xpos.x <- xpos[rk.x]
  xpos.y <- xpos[rk.y]
  ypos.x <- ypos[rk.x]
  ypos.y <- ypos[rk.y]

  ## call routine 
  z <- .C("cor_spatial",
          xpos.x = as.double(xpos.x),
          xpos.y = as.double(xpos.y),
          ypos.x = as.double(ypos.x),
          ypos.y = as.double(ypos.y),
          bars = as.double(bars),
          xpos = as.double(xpos),
          ypos = as.double(ypos),
          dims = as.integer(dims),
          cor = as.double(0),
          var = as.double(0))
  
  ## creating output object
  x <- z$cor
  attr(x, "variance") <- z$var
  x
}
