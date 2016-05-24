convex.hull<-function(tri.obj,plot.it=FALSE, add=FALSE,...)
{
  if(!inherits(tri.obj,"tri"))
    stop("tri.obj must be of class \"tri\"")
  tborder<-c(rep(-1,tri.obj$n))
  storage.mode(tborder)<-"integer"
  ans<-.Fortran("bnodes",
                 as.integer(tri.obj$n),
                 as.integer(tri.obj$tlist),
                 as.integer(tri.obj$tlptr),
                 as.integer(tri.obj$tlend),
                 tborder=as.integer(tborder),
                 nb=as.integer(0),
                 na=as.integer(0),
                 nt=as.integer(0),
                 PACKAGE = "tripack")
  ret<-list(x=tri.obj$x[ans$tborder[ans$tborder>0]],
            y=tri.obj$y[ans$tborder[ans$tborder>0]],
            i=seq(1,tri.obj$n)[ans$tborder[ans$tborder>0]])
  if(plot.it)
    {
      if (!add)
        {
          plot.new()
          plot.window(range(ret$x), range(ret$y), "")
        }
      lines(cbind(ret$x,ret$x[1]),cbind(ret$y,ret$y[1]), ...)
      invisible(ret)
    }
  else
    ret
}
