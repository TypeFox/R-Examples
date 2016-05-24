#  Part of the statnet package, http://statnetproject.org
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) in
#    http://statnetproject.org/attribution
#
#  Copyright 2014 the statnet development team
######################################################################

# utility methods for tsna stuff

# create a tree network from the results of a path search
as.network.tPath<-function(x,...){
  if(!is.tPath(x)){
    stop("as.network.tPath requires an object of class 'tPath'")
  }
  distance<-x$tdist
  previous<-x$previous
  tree<-network.initialize(length(distance),directed=TRUE)
  vids<-which(distance<Inf)
  for(v in seq_along(vids)){
    if(previous[vids[v]]!=0){ # source vertex will have previous id of 0, so out of range
      fromId<-previous[vids[v]]
      add.edges.active(tree,tail=fromId,head=vids[v],onset=distance[vids[v]],terminus=Inf)
    }
  }
  tree%v%'tdist'<-distance
  tree%v%'gsteps'<-x$gsteps
  return(tree)
}

is.tPath<-function(x){
  if('tPath'%in%class(x)){
    return(TRUE)
  } else {
    return(FALSE)
  }
}

plot.tPath<-function(x,edge.col="red",
                     edge.label.col=edge.col,
                     edge.lwd=10,
                     edge.label.cex=0.7,
                     displaylabels=TRUE,
                     displayisolates=FALSE,
                     jitter=FALSE,
                     vertex.lwd=(x$gsteps==0)*4+1,
                     vertex.cex=(x$gsteps==0)*1.5,
                     vertex.col=NA,...){
  tree<-as.network(x)
  edgeTimes<-sapply(get.edge.activity(tree),'[',1)
  # plot the tree as an overlay
  plot.network(tree,
               displaylabels=displaylabels,
               displayisolates=displayisolates,
               edge.lwd=edge.lwd,
               edge.col=edge.col,
               edge.label=edgeTimes,
               edge.label.col=edge.label.col,
               edge.label.cex=edge.label.cex,
               vertex.lwd=vertex.lwd,
               vertex.cex=vertex.cex,
               vertex.border=edge.col,
               vertex.col=vertex.col,
               jitter=jitter,...)
 
}

# plot a network with a hilited path
# and some sensible defaults
plotPaths<-function(nd,paths,
                    path.col=rainbow(length(paths),alpha=0.5),
                    displaylabels=TRUE, coord = NULL,...){
  # plot the network normally and save coords (if not already passed in)
  coords<-plot.network(nd,displaylabels=displaylabels,coord=coord,...)
  # check if it is a single path or a list of paths
  if (is.tPath(paths)){
    plot(paths,
         coord=coords,
         edge.col=path.col[1],
         new=FALSE, # to make sure it overplots on existing
         displaylabels=FALSE,
         displayisolates=TRUE, # need to include isolates or it messages up the line scaling becauses sizes to area of vertices actually drawn
         ...)
  } else {
    path.col<-rep(path.col,length(paths))
    for (p in 1:length(paths)){
      # create another network that is the tree and overplot it
      plot(paths[[p]],
           coord=coords,
           edge.col=path.col[p],
           new=FALSE, # make sure it overplotts on existing
           displaylabels=FALSE,
           displayisolates=TRUE, # need to include isolates or it messages up the line scaling becauses sizes to area of vertices actually drawn
           ...)
    }
  }
  invisible(coords)
}

# helper function to determine an appropriate finite start and
# end range for the network using net.obs.period if it exists
get_bounds<-function(nd){
  bounds<-c(0,1)
  obs<-nd%n%'net.obs.period'
  if (!is.null(obs)){
    bounds<-range(unlist(obs$observations))
  } else {
    times<-get.change.times(nd)
    # its possible that network has only INFs
    if(length(times)>0){
      bounds<-range(times)
    }
  }
  return(bounds)
}