## R interface to seedfill ("color" fill with double precision numbers)
## for matrices, useful for grid-based models and spatial statistics

seedfill <- function(z, x=1, y=1, fcol=0, bcol=1, tol=1e-6) {
  if (!is.matrix(z))     stop("z must be a matrix")
  if (!is.numeric(x))    stop("x must be numeric")
  if (!is.numeric(y))    stop("y must be numeric")
  if (!is.numeric(fcol)) stop("fcol must be numeric")
  if (!is.numeric(bcol)) stop("bcol must be numeric")
  if (!is.numeric(tol))  stop("tol must be numeric")
  n <- dim(z)[1]
  m <- dim(z)[2]
  ## the algorithm used has a well-known problem
  ## with cells which have already the fill color
  ## fix: we fill with a color not contained within z
  ffcol <- .Machine$double.xmax
  z <-.C("seedfill",
         as.integer(n),
         as.integer(m),
         as.integer(x-1),
         as.integer(y-1),
         z=as.double(z),
         as.double(ffcol),
         as.double(bcol),
         as.double(tol),
         PACKAGE="simecol")$z
  ## fix: replace back ffcol with fcol
  z<-ifelse(z == ffcol, fcol, z)
  return(matrix(z, nrow=n, ncol=m))
}
