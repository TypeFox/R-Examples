triangles<-function(tri.obj)
{
  if(!inherits(tri.obj,"tri"))
    stop("tri.obj must be of class \"tri\"")
  nt<-summary(tri.obj)$nt;
  ans<-.Fortran("trlist",
                 as.integer(tri.obj$nc),
                 as.integer(tri.obj$lc),
                 as.integer(tri.obj$n),
                 as.integer(tri.obj$tlist),
                 as.integer(tri.obj$tlptr),
                 as.integer(tri.obj$tlend),
                 as.integer(9),
                 as.integer(nt),
                 tltri=integer(9*nt),
                 lct=integer(tri.obj$nc),
                 ier=as.integer(0),
                 PACKAGE = "tripack")
  ret<-matrix(ans$tltri,nt,9,byrow=TRUE)
  colnames(ret)<-c("node1","node2","node3","tr1","tr2","tr3","arc1","arc2","arc3")
  ret
}
