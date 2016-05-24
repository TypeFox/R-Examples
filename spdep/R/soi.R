# Copyright 2001 by Nicholas Lewin-Koh
#


soi.graph <- function(tri.nb, coords){
  x <- coords
  if (!is.matrix(x)) stop("Data not in matrix form")
  if (any(is.na(x))) stop("Data cannot include NAs")
  np<-length(tri.nb)
  noedges<-0
  rad<-nearneigh<-rep(0,np)
  neigh<-unlist(tri.nb)  
  noneigh<-unlist(lapply(tri.nb,length))
  g1<-g2<-rep(0,sum(noneigh))
  storage.mode(x) <- "double"
  answ<-.C("compute_soi", np=as.integer(np), from=as.integer(g1),
     to=as.integer(g2), nedges=as.integer(noedges),
     notri.nb=as.integer(noneigh), tri.nb=as.integer(neigh),
     nn=as.integer(nearneigh), 
     circles=as.double(rad), x=x[,1], y=x[,2],
     PACKAGE="spdep")
  answ$from<-answ$from[1:answ$nedges]
  answ$to<-answ$to[1:answ$nedges]
  answ<-list(np=answ$np,nedges=answ$nedges,
             from=answ$from,to=answ$to,circles=answ$circ)
  attr(answ, "call") <- match.call()
  class(answ)<-c("Graph","SOI")
  answ
}
