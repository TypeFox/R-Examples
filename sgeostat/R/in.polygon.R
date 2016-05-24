"in.polygon" <-
  function (x0,y0,x,y)
{
   n <- length(x)
  n0 <- length(x0)
  if(length(y0)!=n0)
    stop("length of x0 and y0 differ!")
  ret <- .Fortran("inplg",
                  as.double(c(x,0)),
                  as.double(c(y,0)),
                  as.integer(n),
                  as.double(x0),
                  as.double(y0),
                  as.integer(n0),
                  inhull = integer(n0),
                  PACKAGE="sgeostat")
  as.logical(ret$inhull)
}
