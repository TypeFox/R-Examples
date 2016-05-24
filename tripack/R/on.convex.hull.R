on.convex.hull<-function(tri.obj,x,y)
{
  if(!inherits(tri.obj,"tri"))
    stop("tri.obj must be of class \"tri\"")
  if(length(x)!=length(y))
    stop("x and y must be of same length")
  n<-length(x)
  if(n==0)
    stop("length of x (resp. y) is 0")
  ans<-.Fortran("onhull",
                as.double(x),
                as.double(y),
                as.integer(n),
                as.double(tri.obj$x),
                as.double(tri.obj$y),
                as.integer(tri.obj$n),
                as.integer(tri.obj$tlist),
                as.integer(tri.obj$tlptr),
                as.integer(tri.obj$tlend),
                onhull=logical(n),
                eps=as.double(1E-15),
                PACKAGE = "tripack")
  ans$onhull
}
