#  File degreenet/R/simyule.R
#  Part of the statnet package, http://statnet.org
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) in
#    http://statnet.org/attribution
#
# Copyright 2003 Mark S. Handcock, University of California-Los Angeles
# Copyright 2007 The statnet Development Team
######################################################################
#
# Find a graph from a Yule distribution
# but not a random graph
#
ryule <- function(n=20,rho=2.5, maxdeg=n-1,maxit=10,verbose=FALSE){
 if (!requireNamespace("igraph", quietly = TRUE) | !requireNamespace("network", quietly = TRUE)) {
  stop('The simyule function requires both the "igraph" and "network" packages to be available.')
 }
 mdeg <- n 
 deg <- 1
 while(mdeg>n-1 | 2*floor(sum(deg)/2) != sum(deg) ){
  if(verbose & sum(deg)>1){
   print(paste("odd number of links =",sum(deg),", so resampling..."))}
  deg <- sample(x=1:maxdeg,size=n,replace=TRUE,prob=dyule(v=rho,x=1:maxdeg))
  mdeg <- max(deg)
  if(verbose & mdeg>n-1){print(paste("max.deg =",mdeg,"> n - 1, so resampling..."))}
  if(verbose){print(table(deg))}
  deg[deg>n-1] <- n-1
  mdeg <- max(deg)
 }
 deg <- deg[order(-deg)]
 sm <- .catchToList(igraph::get.edgelist(igraph::degree.sequence.game(deg,method="vl")))
 iter <- 0
 if(!is.null(sm$error)){
   sm <- .catchToList(igraph::get.edgelist(igraph::degree.sequence.game(deg,method="simple.no.multiple")))
   while(!is.null(sm$error) & iter < maxit){
    while(mdeg>n-1 | 2*floor(sum(deg)/2) != sum(deg) ){
     if(verbose & sum(deg)>1){
      print(paste("odd number of links =",sum(deg),", so resampling..."))}
     deg <- sample(x=1:maxdeg,size=n,replace=TRUE,prob=dyule(v=rho,x=1:maxdeg))
     mdeg <- max(deg)
     if(verbose & mdeg>n-1){print(paste("max.deg =",mdeg,"> n - 1, so resampling..."))}
     if(verbose){print(table(deg))}
     deg[deg>n-1] <- n-1
     mdeg <- max(deg)
    }
    deg <- deg[order(-deg)]
    sm <- .catchToList(igraph::get.edgelist(igraph::degree.sequence.game(deg,method="simple.no.multiple")))
    iter <- iter + 1
   }
 }
 if(iter >= maxit){
  stop('The simyule function failed to form a valid network from the passed degree sequence.')
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
