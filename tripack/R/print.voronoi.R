print.voronoi<-function(x,...)
{
  if(!inherits(x,"voronoi"))
    stop("x must be of class \"voronoi\"")
  cat("voronoi mosaic:\n")
  cat("nodes: (x,y): neighbours (<0: dummy node)\n")
  for (i in 1:length(x$x))
    {
      if(x$node[i]){
        cat(i,": (",x$x[i],",",x$y[i],")",sep="")
        cat(":",x$n1[i],x$n2[i],x$n3[i],"\n",sep=" ")
      }
    }
  cat("dummy nodes: (x,y)\n")
  for (i in 1:length(x$dummy.x))
    {
      cat(i,": (",x$dummy.x[i],",",x$dummy.y[i],")\n",sep="")
    }

}
