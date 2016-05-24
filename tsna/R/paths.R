#  Part of the statnet package, http://statnetproject.org
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) in
#    http://statnetproject.org/attribution
#
#  Copyright 2014 the statnet development team
######################################################################

# functions for evaluating temporal paths in networks

# this is a wrapper function to check args and call the appropriate paths method
tPath<-function(nd,v, 
                 direction=c('fwd','bkwd'),
                 type=c('earliest.arrive', 'latest.depart'),
                 start,end,active.default=TRUE,
                 graph.step.time=0){
  
  if (!is.networkDynamic(nd)){
    stop('to be able to calculate temporal paths, the first argument must be a networkDynamic object')
  }
  if (missing(v) || !is.numeric(v)){
    stop("a 'v' argument with valid vertex ids was not given to specify starting vertex")
  }
  
  if (network.size(nd)>0 && !v%in%seq_len(network.size(nd))){
    stop(" the 'v' argument must be a vertex id in the appropriate range for the network size")
  }
  
  if (!missing(start)&!missing(end)){
    if (!is.null(start)&!is.null(end)){
      if(start>end){
        stop("the time value for the 'start' parameter may not be greater than the 'end' parameter")
      }
    }
  }
  
  direction<-match.arg(direction)
  type<-match.arg(type)
  
  # check and determine starting and ending times
  
  # determine which version to call
  values<-NULL
  if (direction=='fwd'){
    if (type=='earliest.arrive'){
      values <- paths.fwd.earliest(nd=nd,v=v,start=start,end=end,active.default=active.default,graph.step.time=graph.step.time)
    } else if (type=='fewest.steps'){
      values <- paths.fwd.fewest.steps(nd=nd,v=v,start=start,end=end,active.default=active.default,graph.step.time=graph.step.time)
    }
    
  } else {  # direction is backwards
    if (type=='latest.depart'){
      values <- paths.bkwd.latest(nd=nd,v=v,start=start,end=end,active.default=active.default,graph.step.time=graph.step.time)
    } 
  }
  if (is.null(values)){
    stop('unable to calculate ',direction, ' ',type, ' paths because the method is not yet implemented')
  }
  # record meta information about path type
  values$'direction'<-direction
  values$'type'<-type
  # prepent the 'tPath' class
  class(values)<-c('tPath',class(values))
  return(values)
}


# compute the forward-shortest path with a Dijkstra-style earch
# search stats from vertex v
# temporal search is bounded by 'start' and 'end' times. 

# TODO: add option to make direct of edge evaluation explicit?

# this finds an earliest-ending path
paths.fwd.earliest<-function(nd,v,start,end,active.default=TRUE,graph.step.time=0,alter){
  
  # check for negatively valued times to catch error condition in #1135
  changes<-get.change.times(nd)
  if(any(changes<0)){
    warning("paths.fwd.earliest may give correct results if any of the edges have spells with time values less than zero.")
  }
  
  if (missing(start) || is.null(start)){
    # TODO: use obs.period if it exists
    #changes<-get.change.times(nd)
    if(length(changes)>0){
      start<-min(changes)
      # message("'start' parameter was not specified, using value first network change '",start)
    } else {
      # can't use inf, because all distances will be inf
      start<-0
      message("'start' time parameter for paths was not specified, no network changes found,  using start=",start)
    }
  }
  if (missing(end) || is.null(end)){
    # TODO: use obs.period if it exists
    end<-Inf
  }
  
  if (graph.step.time<0){
    stop("'graph.step.time' paramter must be a positive value")
  }
  
  # TODO: self-loop behavior?
  # TODO: multiplex behavior?
  
  dist<-rep(Inf,network.size(nd))     # default all vertices to un-reachable (Inf)
  previous<-rep(0,network.size(nd)) # array used for reconstructing the path
  dist[v]<-0                          # set distance to self to 0
  toCheck<-rep(TRUE,network.size(nd)) # make a list of unchecked vertices
  # begin depth first search loop
  while(sum(toCheck)>0){
    minToCheck<-which.min(dist[toCheck])  # select 'closest' vertex to check
    # have to translate index found back to network index
    u<-which(toCheck)[minToCheck]
    toCheck[u]<-FALSE
    if (dist[u]>= end-start){  #NOTICE:  DISTANCE IS NOT ABSOLUTE TIME  should be end-start
      break;  # no more vertices are reachable from v within time range
    }
    
    nghE<-get.edgeIDs(nd,v=u,neighborhood='out') # check neighbors of u
    for (e in nghE){
      # get vertex index of u's neighbor
      w <- ifelse(nd$mel[[e]]$inl==u,nd$mel[[e]]$outl,nd$mel[[e]]$inl)   
      # we ignore graph hop time
      # so "distance" is how long we have to wait from 'now' until onset of edge
      spls<-nd$mel[[e]]$atl$active
      if (is.null(spls)){ # handle possibly missing activity value, assume always active
        if (active.default){
          dist_u_w<-0+graph.step.time 
        } else {
          dist_u_w<-Inf
        }
        
      } else { # since edge activities are defined ..
        # find the index of the active spell
        splIndex<-spells.hit(needle=c(start+dist[u],end),haystack=spls)
        # if we are using graph.step.time > 0, may need to search for a later spells
        while (splIndex>0){
          # check remaining duration of edge spell is not too short transmission to occur
          if (max(0,(spls[splIndex,1]-start)-dist[u])+graph.step.time > spls[splIndex,2]){
            # search forward to see if there are any more spells active after selected spell
            if(splIndex<nrow(spls)){
              s <- splIndex+1
              splIndex <- -1
              for (s in s:nrow(spls)) {
                if (spells.overlap(c(start+dist[u],end), spls[s, ])) {
                  splIndex<-s
                  break()
                }
              }
            } else {
              splIndex<- -1 # we were already checking the last spell
            }
          } else {
            break() # this spell duration is ok, so keep on going with it
          }
        }
        if (splIndex<0){ # no active spells found so
          dist_u_w<-Inf  # vertex is never reachable in the future / within time bound
        
        } else {
          # otherwise additional distance is the later of 0 or the difference between the 
          # 'current' time and the onset of the edge
          # if we are counting graph steps as part of the distance, distance can't be less than graph step
          dist_u_w<-max(0,(spls[splIndex,1]-start)-dist[u])+graph.step.time #
        }
      }
      dist_v_w <-dist[u]+dist_u_w 
      if (dist_v_w < dist[w]){ # if this new value is shorter, update
        dist[w]<-dist_v_w
        previous[w]<-u
      }
    }
    # check if we have found the target alter, if so, stop searching
    if (!missing(alter)){
      if(all(toCheck[alter]==FALSE)){ # if all alters are checked
        break;
      }
    }
  }
  
  # construct the vector of geodeisc graph steps from previous
  gsteps<-rep(Inf,length(previous))
  curDist<-0
  preV<-v
  # loop down the tree path and mark the distances
  # until all reachable vertices are updated
  while(length(preV)>0){
    gsteps[preV]<-curDist
    preV<-unlist(sapply(preV,function(v){
      which(previous==v)
    }))
    curDist<-curDist+1
  }

  return(list(tdist=dist,previous=previous,gsteps=gsteps, start=start, end=end))
}

# this finds the forward path with the fewest number of intermediate vertices
# THIS IS WRONG, IT RETURNS THE  SHORTEST PATH THAT IS ALSO EARLISET (but may miss shortest paths arriving later) or something like that
paths.fwd.fewest.steps<-function(nd,v,start,end,active.default=TRUE,graph.step.time=0){
  
  if (missing(start) || is.null(start)){
    # TODO: use obs.period if it exists
    changes<-get.change.times(nd)
    if(length(changes)>0){
      start<-min(changes)
      # message("'start' parameter was not specified, using value first network change '",start)
    } else {
      # can't use inf, because all distances will be inf
      start<-0
      message("'start' time parameter for paths was not specified, no network changes found,  using start=",start)
    }
  }
  if (missing(end) || is.null(end)){
    # TODO: use obs.period if it exists
    end<-Inf
  }
  
  if (graph.step.time<0){
    stop("'graph.step.time' paramter must be a positive value")
  }
  
  # TODO: self-loop behavior?
  # TODO: multiplex behavior?
  
  times<-rep(Inf,network.size(nd))     # default all vertices to un-reachable (Inf)
  gsteps<-rep(Inf,network.size(nd))
  previous<-rep(0,network.size(nd)) # array used for reconstructing the path
  times[v]<-0                          # set distance to self to 0
  gsteps[v]<-0
  toCheck<-rep(TRUE,network.size(nd)) # make a list of unchecked vertices
  # begin depth first search loop
  while(sum(toCheck)>0){
    minToCheck<-which.min(gsteps[toCheck])  # select 'closest' vertex to check in terms of gsteps
    # have to translate index found back to network index
    u<-which(toCheck)[minToCheck]
    toCheck[u]<-FALSE
    if (times[u]>= end-start){  #NOTICE:  DISTANCE IS NOT ABSOLUTE TIME  should be end-start
      break;  # no more vertices are reachable from v within time range
    }
    
    nghE<-get.edgeIDs(nd,v=u,neighborhood='out') # check neighbors of u
    for (e in nghE){
      # get vertex index of u's neighbor
      w <- ifelse(nd$mel[[e]]$inl==u,nd$mel[[e]]$outl,nd$mel[[e]]$inl)   
      # we ignore graph hop time
      # so "distance" is how long we have to wait from 'now' until onset of edge
      spls<-nd$mel[[e]]$atl$active
      if (is.null(spls)){ # handle possibly missing activity value, assume always active
        if (active.default){
          time_u_w<-0+graph.step.time
          gsteps_u_w<-1 
        } else {
          time_u_w<-Inf
          gsteps_u_w<-Inf
        }
        
      } else { # since edge activities are defined ..
        # find the index of the active spell
        splIndex<-spells.hit(needle=c(start+times[u],end),haystack=spls)
        # if we are using graph.step.time > 0, may need to search for a later spells
        while (splIndex>0){
          #check remaining duration of edge spell is long enough for transmission to occur
          if (max(0,(spls[splIndex,1]-start)-times[u])+graph.step.time > spls[splIndex,2]){
            #query again to see if there are any spells active after selected spell
            splIndex<-spells.hit(needle=c(spls[splIndex,2],end),haystack=spls)
          } else {
            break() # this spell duration is ok, so keep on going with it
          }
        }
        if (splIndex<0){ # no active spells found so
          gsteps_u_w<-Inf
          time_u_w<-Inf  # vertex is never reachable in the future / within time bound
          
        } else {
          # otherwise additional distance is the later of 0 or the difference between the 
          # 'current' time and the onset of the edge
          # if we are counting graph steps as part of the distance, distance can't be less than graph step
          time_u_w<-max(0,(spls[splIndex,1]-start)-times[u])+graph.step.time #
          gsteps_u_w <- 1 # graph step distances (unweighted)
        }
      }
      gsteps_v_w <-gsteps[u]+gsteps_u_w 
      time_v_w <-times[u]+time_u_w
      if (gsteps_v_w < gsteps[w]){ # if this new GEODESIC DISTANCE is shorter, update distance and times
        gsteps[w]<-gsteps_v_w
        times[w]<-time_v_w
        previous[w]<-u
      }
    }
  }
  
  
  return(list(tdist=times,previous=previous,gsteps=gsteps,start=start,end=end))
}

# compute reverse paths 
# need to start at the end and minimize backwards

paths.bkwd.latest<-function(nd,v,start,end,active.default=TRUE,graph.step.time=0,alter){

  if (missing(end) || is.null(end)){
    # TODO: use obs.period if it exists
    changes<-get.change.times(nd)
    if(length(changes)>0){
      end<-max(changes)
    } else {
      # this means that edges are either always active or never active
      # if they are never active, it won't matter, because all distances other than v will be Inf
      # if they are always active, the longest possible path would be equal to the size of the network
      # times the graph.step.time
      end<-network.size(nd)*graph.step.time
      warning("'end' time parameter for paths was not specified, no network changes found, using 'end' value of ",end)
      
    }
  }
  if (missing(start) || is.null(start)){
    # TODO: use obs.period if it exists
    start<- -Inf
  }
  
  if (graph.step.time<0){
    stop("'graph.step.time' paramter must be a positive value")
  }
  
  # TODO: self-loop behavior?
  # TODO: multiplex behavior?
  
  dist<-rep(Inf,network.size(nd))
  previous<-numeric(network.size(nd)) # array used for reconstructing the path
  dist[v]<-0
  toCheck<-rep(TRUE,network.size(nd))
  while(sum(toCheck)>0){
    minToCheck<-which.min(dist[toCheck])
    # have to translate index found back
    u<-which(toCheck)[minToCheck]
    toCheck[u]<-FALSE
    if (dist[u]>= end-start){
      break;  # no more vertices are reachable from v within time range
    }
    
    # we are going backwards, so use 'in' edges instead of 'out'
    nghE<-get.edgeIDs(nd,v=u,neighborhood='in') # check neighbors of u
    for (e in nghE){
      w <- ifelse(nd$mel[[e]]$inl==u,nd$mel[[e]]$outl,nd$mel[[e]]$inl)
      # we ignore graph hop time
      # so "distance" is how long we have to wait from 'now' until onset of edge
      spls<-nd$mel[[e]]$atl$active
      if (is.null(spls)){ # handle possibly missing activity value, assume always active
        if (active.default){
          dist_u_w<-0+graph.step.time #TODO: need to check active default here to know if returning Inf or dist[w]
        } else {
          dist_u_w<- Inf
        }
        
      } else {
        # can't use spells.hit because it returns earliest spell, not latest
        splIndex<- -1
        # loop backwards over spells so we find latest first
        for (s in nrow(spls):1) {
          if (spells.overlap(c(start,end-dist[u]), spls[s, ])) {
            splIndex<-s
            break
          }
        }
        
        # if we are using graph.step.time > 0, may need to search for a later spells
        while (splIndex>0){
          #check remaining duration of edge spell is long enough for transmission to occur
                     #now     # start of spell
          if (  ((end-dist[u]) - spls[splIndex,1]) < graph.step.time ){
            # spell was not long enough so query again to see if there were any other 
            # spells active earlier than the one we just found
            # (can't use spells.hit because it returns earliest spell, not latest)
            splIndex<-splIndex-1
            # loop backwards over spells so we find latest first
            while (splIndex > 0) {
              if (spells.overlap(c(start,end-dist[u]), spls[splIndex, ])) {
                break
              }
              splIndex<-splIndex-1
            }
          } else {
            break() # this spell duration is ok, so keep on going with it
          }
        }
        
        if (splIndex>0){
          # distance is the later of dist[u] or the terminus of the edge,
          # but distance can't be less than grap.step.time
          dist_u_w<-max(0,(spls[splIndex,2]-end)*-1-dist[u])+graph.step.time
        } else {
          # otherwise vertex is never reachable in the future / within time bound
          dist_u_w<- Inf  
        }
      }
      dist_v_w <-dist[u]+dist_u_w 
      if (dist_v_w < dist[w]){ # if this new value is shorter, update
        dist[w]<-dist_v_w
        previous[w]<-u
      }
    } # end edge neighbor compare loop
    # check if we have found the target alter, if so, stop searching
    if (!missing(alter)){
      if(all(toCheck[alter]==FALSE)){ # if all alters are checked
        break;
      }
    }
  } # end search loop
  
  # construct the vector of geodeisc graph steps from previous
  gsteps<-rep(Inf,length(previous))
  curDist<-0
  preV<-v
  # loop down the tree path and mark the distances
  # until all reachable vertices are updated
  while(length(preV)>0){
    gsteps[preV]<-curDist
    preV<-unlist(sapply(preV,function(v){
      which(previous==v)
    }))
    curDist<-curDist+1
  }
  # TODO: we are measuring distance backwards from the end
  # so need to flip distance measure
  
  return(list(tdist=dist,previous=previous, gsteps=gsteps,start=start,end=end))
}

# this version tries to minimize the distance of the latest time forward
# calculate the latest ending
# THIS DOES NOT WORK (and I think it can't, due to longest path problem)
paths.fwd.latestBAD<-function(nd,v,start,end,active.default=TRUE,graph.step.time=0){
  
  if (!is.networkDynamic(nd)){
    stop('to be able to calculate forward paths, the first argument must be a networkDynamic object')
  }
  if (missing(v) || !is.numeric(v)){
    stop("a 'v' argument with valid vertex ids was not given to specify starting vertex")
  }
  if (missing(end)){
    # TODO: use obs.period if it exists
    changes<-get.change.times(nd)
    if(length(changes)>0){
      end<-max(changes)
      # message("'start' parameter was not specified, using value first network change '",start)
    } else {
      stop("'end' time parameter for paths was not specified, no network changes found")
    }
  }
  if (missing(start)){
    # TODO: use obs.period if it exists
    start<- -Inf
  }
  
  
  # TODO: self-loop behavior?
  # TODO: multiplex behavior?
  
  dist<-rep(Inf,network.size(nd))
  previous<-numeric(network.size(nd)) # array used for reconstructing the path
  dist[v]<-0
  toCheck<-rep(TRUE,network.size(nd))
  while(sum(toCheck)>0){
    minToCheck<-which.min(dist[toCheck])
    # have to translate index found back
    u<-which(toCheck)[minToCheck]
    toCheck[u]<-FALSE
    if (dist[u]>= end-start){
      break;  # no more vertices are reachable from v within time range
    }
    
    # we are going forwards, so use 'out' edges instead of 'in'
    nghE<-get.edgeIDs(nd,v=u,neighborhood='out') # check neighbors of u
    for (e in nghE){
      w <- ifelse(nd$mel[[e]]$inl==u,nd$mel[[e]]$outl,nd$mel[[e]]$inl)   
      # we ignore graph hop time
      # so "distance" is how long we have to wait from 'now' until onset of edge
      spls<-nd$mel[[e]]$atl$active
      if (is.null(spls)){ # handle possibly missing activity value, assume always active
        if (active.default){
          dist_u_w<-0 #TODO: need to check active default here to know if returning Inf or dist[w]
        } else {
          dist_u_w<- Inf
        }
        
      } else {
        # can't use spells.hit because it returns earliest spell, not latest
        splIndex<- -1
        # loop backwards over spells so we find latest first
        for (s in nrow(spls):1) {
          if (spells.overlap(c(start+dist[u],end), spls[s, ])) {
            splIndex<-s
            break
          }
        }
        
        if (splIndex<0){
          dist_u_w<- Inf  # vertex is never reachable in the future / within time bound
        } else {
          # otherwise distance is the later of dist[u] or the terminus of the edge
          dist_u_w<-max(0,(spls[splIndex,2]-end)*-1-dist[u]+graph.step.time)
        }
      }
      dist_v_w <-dist[u]+dist_u_w 
      if (dist_v_w < dist[w]){ # if this new value is shorter, update
        dist[w]<-dist_v_w
        previous[w]<-u
      }
    }
  }
  
  # TODO: we are measuring distance backwards from the end
  # so need to flip distance measure
  
  return(list(tdist=dist,previous=previous))
}


# this actually calculates the forward path to find the reachable set, 
# then calculates the latest *backward path* from each reached vertex, 
# and then reverse the values
paths.fwd.latest<-function(nd,v,start,end,active.default=TRUE,graph.step.time=0)
{
  warning("paths.fwd.latest has not be fully tested and may not be correct")
  if (missing(start)){
    # TODO: use obs.period if it exists
    changes<-get.change.times(nd)
    if(length(changes)>0){
      start<-min(changes)
      # message("'start' parameter was not specified, using value first network change '",start)
    } else {
      # can't use inf, because all distances will be inf
      start<-0
      message("'start' time parameter for paths was not specified, no network changes found,  using start=",start)
    }
  }
  if (missing(end)){
    # TODO: use obs.period if it exists
    changes<-get.change.times(nd)
    if(length(changes)>0){
      end<-max(changes)
      # message("'start' parameter was not specified, using value first network change '",start)
    } else {
      # can't use inf, because all distances will be inf
      end<-Inf
      message("'end' time parameter for paths was not specified, no network changes found,  using end=",end)
    }
  }
  latest<-rep(Inf, network.size(nd))
  previous<-as.list(rep(0,network.size(nd)))
  # find the earliest path to reachable vertices
  fwdReachable<-which(paths.fwd.earliest(nd=nd,v=v,start=start,end=end,active.default=active.default,graph.step.time=graph.step.time)$tdist<Inf)
  
  # for each reachable vertex, find the latest return path
  
  latestResults<-lapply(fwdReachable,function(w){
    backwards<-paths.bkwd.latest(nd=nd,v=w,start=start,end=end,active.default=active.default,graph.step.time=graph.step.time)
    return(list(backwards$tdist[v],backwards$previous))
  })
  latest[fwdReachable]<-sapply(latestResults,'[[',1)
  previous[fwdReachable]<-lapply(latestResults,'[[',2)
  
  # because distances are from the end, need to reverse by subtracting the end time
  latest[fwdReachable]<-end-latest[fwdReachable]
  return(list(tdist=latest,previous=previous))
}



# testing an approximate / stochastic method for computing arrival probabilities
# works by randomly sampling fwd paths from vertex v
#Currently a number of \code{tries} are made from the starting vertex. For each try, a random hop forward in time #is made, neighbors at that time are discovered, a random choice is made between each neighbor and remaining on #the vertex, a new random hop is made, etc.  Counts are are made of the number of times each vertex is reached, and this value is returned, divided by the number of trials  Due to the stochastic nature of the algorithm, it 
# WILL NOT necessarily find paths to all reachable vertices. 
paths.fwd.approx<-function(nd,v,tries=network.size(nd)*100,mean.hop.dur=1, start,end){
  if (missing(start) | missing(end)){
    bounds<-range(get.change.times(nd,vertex.activity=FALSE,vertex.attribute.activity=FALSE,edge.attribute.activity=FALSE,network.attribute.activity=FALSE))
  }
  if(missing(start)){
    start<-bounds[1]
  }
  if(missing(end)){
    end<-bounds[2]
  }
  if (is.infinite(start)){ 
    stop("Can not evaluate paths over a time interval with an infinite start value")
  }
  if (is.infinite(end)){
    stop("Can not evaluate paths over a time interval with an infinite end value")
  }
  if (is.infinite(mean.hop.dur)){
    stop("Can not evaluate paths over a time interval with an infinite mean.hop.dur")
  }
  
  # TODO: this could be parallelized
  arrivals<-replicate(tries,{
    # start at v at time start
    currentV<-v
    currentTime<-start+rexp(1,1/mean.hop.dur)
    # repeat until run out of time
    while (currentTime<end){
      # find reachable neighbors at that time
      ngs<-get.neighborhood.active(nd,v=currentV,at=currentTime,type='out')
      ngs<-c(currentV,ngs) # include 0 for no-branch option
      # randomly pick a neighbor or self to hop again
      currentV<-sample(ngs,1,replace=TRUE)
      # hop a random amount forward in time 
      #TODO: how to correctly hop randomly
      currentTime<-currentTime+rexp(1,1/mean.hop.dur)
    }
    return(currentV)
  })
  arrivalCounts<-tabulate(arrivals,nbins=network.size(nd))
  return(arrivalCounts/tries)
}

