tri.find<-function(tri.obj,x,y)
{
  if(!inherits(tri.obj,"tri"))
    stop("tri.obj must be of class \"tri\"")
  ans<-.Fortran("trfind",
                 as.integer(1),
                 as.double(x),
                 as.double(y),
                 as.integer(tri.obj$n),
                 as.double(tri.obj$x),
                 as.double(tri.obj$y),
                 as.integer(tri.obj$tlist),
                 as.integer(tri.obj$tlptr),
                 as.integer(tri.obj$tlend),
                 i1=as.integer(0),
                 i2=as.integer(0),
                 i3=as.integer(0),
                 PACKAGE = "tripack")
  list(i1=ans$i1,i2=ans$i2,i3=ans$i3)
}
