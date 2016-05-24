print.tri<-function(x,...)
{
  if(!inherits(x,"tri"))
    stop("x must be of class \"tri\"")
  tnabor<- integer(x$tlnew)
  nnabs <- integer(x$n)
  nptr <- integer(x$n)
  nptr1 <- integer(x$n)
  nbnos <- integer(x$n)
  ans<-.Fortran("troutq",
                 as.integer(x$nc),
                 as.integer(x$lc),
                 as.integer(x$n),
                 as.double(x$x),
                 as.double(x$y),
                 as.integer(x$tlist),
                 as.integer(x$tlptr),
                 as.integer(x$tlend),
                 as.integer(6),
                 nnabs=as.integer(nnabs),
                 nptr=as.integer(nptr),
                 nptr1=as.integer(nptr1),
                 tnabor=as.integer(tnabor),
                 nbnos=as.integer(nbnos),
                 na=as.integer(0),
                 nb=as.integer(0),
                 nt=as.integer(0),
                 PACKAGE = "tripack")
  cat("triangulation nodes with neigbours:\n")
  cat("node: (x,y): neighbours\n")
  for (i in 1:x$n)
    {
      cat(i,": (",x$x[i],",",x$y[i],") [",ans$nnabs[i],"]",sep="")
      cat(":",sort(ans$tnabor[ans$nptr[i]:ans$nptr1[i]]),"\n",sep=" ")
    }
  cat("number of nodes:",x$n,"\n")
  cat("number of arcs:",ans$na,"\n")
  cat("number of boundary nodes:",ans$nb,"\n")
  cat("boundary nodes: ",ans$nbnos[1:ans$nb], "\n", sep=" ")
  cat("number of triangles:",ans$nt,"\n")
  cat("number of constraints:",x$nc,"\n")
}
