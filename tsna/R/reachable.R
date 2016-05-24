#  Part of the statnet package, http://statnetproject.org
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) in
#    http://statnetproject.org/attribution
#
#  Copyright 2014 the statnet development team
######################################################################

# find the set of vertices reachable from seed vertex within time range

forward.reachable<-function(nd,v,start=NULL,end=NULL,per.step.depth=Inf ){
  if(!is.networkDynamic(nd)){
    stop('the first argument to forward.reachble must be a networkDynamic object')
  }
  if(missing(v) || !is.numeric(v)){
    stop('v argument to forward.reachable must be a vector of valid numeric vertex ids')
  }
  if(max(v)>network.size(nd) | min(v)<1){
    stop('v argument to forward.reachable must be a vector of numeric vertex ids within the range of the network size')
  }
  
  # set the interval to be whatever the observed changes are
  times <-get.change.times(nd,vertex.attribute.activity=FALSE,edge.attribute.activity=FALSE,network.attribute.activity=FALSE)

  if(length(times)==0){
    times <-c(0,Inf)
  }
  if(is.null(start)){
    #start<-min(times)
    start<--Inf
  }
  if(is.null(end)){
    #end<-max(times)
    end<-Inf
  }
  # trim times to desired range, making sure to include start and end
  times<-unique(c(start,times[times>=start]))
  times<-unique(c(times[times<=end],end))
  
  distance<-rep(Inf,network.size(nd))
  distance[v]<-times[1]
  
  #TODO: could probably skip all times earlier that the active times in v?
  reached <-v
  for(t in 1:(length(times)-1)){
    # BFS to depth rate
    new<-reached
    # how long until next change?
    duration<-times[t+1]-times[t]
    
    # remove any in the set we've already visited
    if (duration>0){
      d<-1 # we are assuming all geodesic steps count as 1, harder if we calc per edge..
      # keep searching until we reach bounds or run out of verts to find
      # also stop if we find all the vertices
      while(d <= per.step.depth*duration & length(reached)<network.size(nd)){
        ngs<-unlist(unique(sapply(new,function(i){get.neighborhood.active(nd,v=i,at=times[t],type='out')})))
        new<-setdiff(ngs,reached)
        if(length(new)==0){
          break # no more verts to find
        }
        distance[new]<-times[t]
        reached<-c(reached,new)
        d<-d+1
      }
    }
  }
  return(reached)
}

# this was implemented with when.next.edge.change with the idea of reducing the
# search space at each time step, but it is currently slower than the implementation above
# and paths.fwd.earliest is dramatically faster than both
# forward.reachable2<-function(nd,v,start=NULL,end=NULL,interval='changes',per.step.depth=Inf ){
#   if(!is.networkDynamic(nd)){
#     stop('the first argument to forward.reachble must be a networkDynamic object')
#   }
#   if(missing(v) || !is.numeric(v)){
#     stop('v argument to forward.reachable must be a vector of valid numeric vertex ids')
#   }
#   if(max(v)>network.size(nd) | min(v)<1){
#     stop('v argument to forward.reachable must be a vector of numeric vertex ids within the range of the network size')
#   }
#   
#   # if start or end is missing set the interval to be whatever the observed changes are
#   if (is.null(start)){
#     start <-min(get.change.times(nd,vertex.attribute.activity=FALSE,edge.attribute.activity=FALSE,network.attribute.activity=FALSE))
#   }
#   if (is.null(end)){
#     end <-max(get.change.times(nd,vertex.attribute.activity=FALSE,edge.attribute.activity=FALSE,network.attribute.activity=FALSE))
#   }
#     
#   # lets not loop for ever!
#   if (is.infinite(start) | is.infinite(end)){
#     stop("start and end values cannot be infinite because search will not terminate")
#   }
#   
#   #TODO: could probably skip all times earlier that the active times in v?
#   reached <-v
#   now<-start
#   while(now<end){
#     # BFS to depth rate
#     new<-reached
#     # how long until next change?
#     nextTime<-when.next.edge.change(nd,at=now,v=reached)
#     duration<-nextTime-now
#     
#     # remove any in the set we've already visited
#     if (duration>0){
#       d<-1 # we are assuming all geodesic steps count as 1, harder if we calc per edge..
#       # keep searching until we reach bounds or run out of verts to find
#       # also stop if we find all the vertices
#       while(d <= per.step.depth*duration & length(reached)<network.size(nd)){
#         ngs<-unlist(unique(sapply(new,function(i){get.neighborhood.active(nd,v=i,at=now,type='out')})))
#         new<-setdiff(ngs,reached)
#         if(length(new)==0){
#           break # no more verts to find
#         }
#         reached<-c(reached,new)
#         d<-d+1
#       }
#     }
#     # update time
#     now<-nextTime
#   }
#   return(reached)
# }



# compute the sets of vertices reachable from each vertex on the graph
tReach<-function(nd,direction=c('fwd','bkwd'),sample=network.size(nd), 
                 start, end, graph.step.time=0){
  
  if (!is.networkDynamic(nd)){
    stop("the first argument must be a networkDynamic object in order to calculate reachable set sizes")
  }
  
  # if sample equals network size, choose seeds so they will be in order
  if(sample==network.size(nd)){
    seeds<-seq_len(network.size(nd))
  } else {
    # otherwise, sample (with replacement)
    seeds<-sample.int(network.size(nd),size=sample)
  }
  
  if (missing(start)) {
    start=NULL
  }
  if (missing(end)) {
    end <- NULL
  }
  
  direction<-match.arg(direction)
  type<-'earliest.arrive'
  if (direction=='bkwd'){
    type<-'latest.depart'
  }
  sizes<-sapply(seeds,function(v){
    sum(tPath(nd,v=v,direction=direction,type=type,start=start,end=end,graph.step.time=graph.step.time)$tdist<Inf)
    })
  return(sizes)
}



