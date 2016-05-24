cells<-function(voronoi.obj)
{
  if(!inherits(voronoi.obj,"voronoi"))
    stop("voronoi.obj must be of class \"voronoi\"")

  tri <- voronoi.obj$tri
  
  tnabor<- integer(tri$tlnew)
  nnabs <- integer(tri$n)
  nptr <- integer(tri$n)
  nptr1 <- integer(tri$n)
  nbnos <- integer(tri$n)
  ans<-.Fortran("troutq",
                 as.integer(tri$nc),
                 as.integer(tri$lc),
                 as.integer(tri$n),
                 as.double(tri$x),
                 as.double(tri$y),
                 as.integer(tri$tlist),
                 as.integer(tri$tlptr),
                 as.integer(tri$tlend),
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

  ret <- NULL
  for(i in 1:tri$n){
    vs <- voronoi.findvertices(i, voronoi.obj)
    if(length(vs)>0){
      center <- c(tri$x[i],tri$y[i])
      neighbours <- sort(ans$tnabor[ans$nptr[i]:ans$nptr1[i]])
      nodes <- rbind(voronoi.obj$x[vs],voronoi.obj$y[vs])
      rownames(nodes) <- c("x","y")
      area <- voronoi.polyarea( voronoi.obj$x[vs], voronoi.obj$y[vs])
      ret[[i]] <- list(cell=i,center=center,
                       neighbours=neighbours,
                       nodes=nodes,area=area)
    } else {
      center <- c(tri$x[i],tri$y[i])
      neighbours <- sort(ans$tnabor[ans$nptr[i]:ans$nptr1[i]])
      nodes <- NA # should better return at least the non-dummy nodes
      area <- NA
      ret[[i]] <- list(cell=i,center=center,
                       neighbours=neighbours,
                       nodes=nodes,area=area)
    }
  }
  
  ret
}

