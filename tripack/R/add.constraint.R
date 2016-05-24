add.constraint<-function(tri.obj,cstx,csty,reverse=FALSE)
{
  if(!inherits(tri.obj,"tri"))
    stop("tri.obj must be of class \"tri\"")
  nt<-summary(tri.obj)$nt;
  lcst<-length(cstx)
  if(reverse)
    {
      cstx<-cstx[lcst:1]
      csty<-csty[lcst:1]
    }
  if(length(csty)!=lcst)
    stop("length of cstx and csty differ")
  # don't modify tri.obj
  tri.obj1<-tri.obj
  if(tri.obj1$nc==0)
    {
      tri.obj1$nc<-tri.obj1$nc+1
      tri.obj1$lc[1]<-tri.obj1$n+1
      tri.obj1$lc[2]<-tri.obj1$n+lcst+1
    }
  else
    {
      tri.obj1$nc<-tri.obj1$nc+1
      tri.obj1$lc[tri.obj1$nc]<-tri.obj1$n+1
      tri.obj1$lc[tri.obj1$nc+1]<-tri.obj1$n+lcst+1
    }
  n1<-tri.obj1$n+1
  n2<-tri.obj1$n+lcst
  tri.obj1$x[n1:n2]<-cstx
  tri.obj1$y[n1:n2]<-csty
  tri.obj1$n<-tri.obj1$n+lcst
  # generate a triangulation with the additional nodes:
  # (we need an updated tlist,tlptr and tlend)
  tri.obj2<-tri.mesh(tri.obj1$x,tri.obj1$y)
  ans<-.Fortran("addcst",
                 as.integer(tri.obj1$nc),
                 as.integer(tri.obj1$lc),
                 as.integer(tri.obj1$n),
                 as.double(tri.obj1$x),
                 as.double(tri.obj1$y),
                 as.integer(2*(tri.obj1$n-3)),
                 integer(2*(tri.obj1$n-3)),
                 tlist=as.integer(tri.obj2$tlist),
                 tlptr=as.integer(tri.obj2$tlptr),
                 tlend=as.integer(tri.obj2$tlend),
                 ier=as.integer(0),
                 PACKAGE = "tripack")
  if(ans$ier==0)
    {
      ret<-list(n=tri.obj1$n,x=tri.obj1$x,y=tri.obj1$y,
                tlist=ans$tlist,tlptr=ans$tlptr,
                tlend=ans$tlend,tlnew=tri.obj2$tlnew,
                nc=tri.obj1$nc,lc=tri.obj1$lc,call=match.call())
    }
  else
    {
      switch(ans$ier,
             stop("nc, n or lc[i] out of range"),
             stop("working array to small"),
             stop("invalid triangulation or collinear nodes on convex hull"),
             stop("intersecting constraint arcs"),
             stop("constraint region contains a node\nmay be you should try \"reverse=TRUE\" to invert the orientation \nof the constraint boundary"),
             )
      stop(paste("error",ans$ier,"in addcst"))
    }
                  
  class(ret)<-"tri"
  invisible(ret)
}
