#  File R/edgelist.R in package network, part of the Statnet suite
#  of packages for network analysis, http://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  http://statnet.org/attribution
#
#  Copyright 2003-2015 Statnet Commons
#######################################################################

# the edgelist functions have been copied in from ergm
as.edgelist <- function(x, ...){
  UseMethod("as.edgelist")
} 


# convert a network into an ergm-style sorted edgelist using
# as.edgelist.matrix and as.matrix.network.edgelist
as.edgelist.network <- function(x, attrname = NULL, as.sna.edgelist = FALSE, inverted = NULL, ...){
  
  # if inverted is given, use it, otherwise check for network property named inverted, otherwise false
  inverted<-ifelse(!is.null(inverted),inverted, ifelse(!is.null(x%n%"inverted"),x%n%"inverted", FALSE))
  as.edgelist(as.matrix.network.edgelist(x, attrname=attrname,
                                         as.sna.edgelist=as.sna.edgelist,...), 
                                         n=network.size(x), 
                                         directed=is.directed(x), 
                                         bipartite=ifelse(is.bipartite(x),x%n%"bipartite",FALSE), 
                                         loops=has.loops(x), 
                                         inverted=inverted,
                                         vnames=network.vertex.names(x))
}



as.edgelist.matrix <- function(x, n, directed=TRUE, bipartite=FALSE, loops=FALSE, inverted=FALSE, vnames=seq_len(n),...){
  if(!directed) x[,1:2] <- cbind(pmin(x[,1],x[,2]),pmax(x[,1],x[,2]))
  if(!loops) x <- x[x[,1]!=x[,2],,drop=FALSE]
  if(bipartite) x <- x[(x[,1]<=bipartite)!=(x[,2]<=bipartite),,drop=FALSE]
  x <- unique(x[order(x[,1],x[,2]),,drop=FALSE])
  attr(x,"n") <- n
  attr(x,"vnames")<- vnames
  attr(x,"directed") <- directed
  attr(x,"bipartite") <- bipartite
  attr(x,"loops") <- loops
  attr(x,"inverted") <- inverted
  class(x)<-c('edgelist',class(x))
  x
}

is.edgelist<-function(x){
  if ('edgelist'%in%class(x)){
    return(TRUE)
  }
  return(FALSE)
}


