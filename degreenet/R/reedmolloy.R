#  File degreenet/R/reedmolloy.R
#  Part of the statnet package, http://statnetproject.org
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) in
#    http://statnetproject.org/attribution
#
# Copyright 2003 Mark S. Handcock, University of California - Los Angeles
# Copyright 2007 The statnet Development Team
######################################################################
#
# Find a graph from a given degree sequence
#
reedmolloy <- function(deg, maxit=10,
                       verbose=TRUE){
# sm is the list of edges (head, tail) counting from 1, 2, ...
 if (!requireNamespace("igraph", quietly = TRUE) | !requireNamespace("network", quietly = TRUE)) {
  stop('The reedmolloy function requires both the "igraph" and "network" packages to be available.')
 }
 sm <- .catchToList(igraph::get.edgelist(igraph::degree.sequence.game(deg,method="vl")))
 iter <- 0
 if(!is.null(sm$error)){
  sm <- .catchToList(igraph::get.edgelist(igraph::degree.sequence.game(deg,method="simple.no.multiple")))
  while(!is.null(sm$error) & iter < maxit){
   jitterdeg <- sample(seq_along(deg),size=2,prob=deg)
   deg[jitterdeg] <- deg[jitterdeg] + 2*(runif(2)>0.5)-1
   deg[deg==0] <- 2
   sm <- .catchToList(igraph::get.edgelist(igraph::degree.sequence.game(deg,method="simple.no.multiple")))
   iter <- iter + 1
  }
 }
 if(iter >= maxit){
   stop('The reedmolloy function failed to form a valid network from the passed degree sequence.')
 }
# pnode <- sample(1:length(deg))
 smn <- network::network.initialize(length(deg), directed=FALSE)
# smn <- add.edges(x=smn,tail=as.list(pnode[sm[,1]+1]),head=as.list(pnode[sm[,2]+1]))
#   tail=as.list(length(deg)-sm[,1]+1),
#   head=as.list(length(deg)-sm[,2]+1))
 smn <- network::add.edges.network(x=smn,
   tail=as.list(sm$value[,1]),
   head=as.list(sm$value[,2]))
 return(smn)
}
