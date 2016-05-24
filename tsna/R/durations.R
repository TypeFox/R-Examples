#  Part of the statnet package, http://statnetproject.org
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) in
#    http://statnetproject.org/attribution
#
#  Copyright 2014 the statnet development team
######################################################################

# functions for calculating and estimating edge durations

edgeDuration<-function(nd,mode=c('duration','counts'),subject=c('edges','spells','dyads'),e=seq_along(nd$mel), start=NULL, end=NULL,active.default=TRUE){
  del<-as.data.frame.networkDynamic(nd,e=e,start=start,end=end,active.default=active.default)
  # if looking per edge, group by edge id
  subject<-match.arg(subject)
  mode<-match.arg(mode)
  
  # if network has no edges, return nothing
  if (nrow(del)==0){
    return(numeric(0))
  }
  
  # determine function type for aggregation
  aggFun<-'sum' # function to use for aggregation sum = 'duration'
  if (mode=='counts'){
    aggFun<-'length'  # function to use to count events
  }
  
  # determine unit of aggragation
  if (subject=='edges'){
    del<-aggregate.data.frame(del[,c('duration','edge.id')],by=list(edges=del$edge.id),FUN=aggFun)
  } else if (subject=='dyads'){
  # if looking per dyad, group by tail,head pair
    if (is.hyper(nd)){
      stop('dyad-based comparison is not appropriate for hypergraphic networks')
    }
    del<-aggregate.data.frame(del[,c('duration','tail','head')],by=list(dyads=del$tail,del$head),FUN=aggFun)
  } else {
  # if looking at spells, just use the raw frame
    del<-aggregate.data.frame(del[,'duration',drop=FALSE],by=list(seq_len(nrow(del))),FUN=aggFun) 
  }
    
  return(del$duration)
  
}


tEdgeFormation<-function(nd, start, end, time.interval=1, result.type=c('count','fraction'),include.censored=FALSE){
    result.type<-match.arg(result.type)
    if(missing(start) | missing(end)){
      times <- get.change.times(nd)
      if (length(times) == 0) {
        warning("network does not appear to have any time range information and start and end parameters not provided. Using start=0 end=0")
        start = 0
        end = 0
      } else {
        times[times == Inf] <- NA
        times[times == -Inf] <- NA
        start = min(times, na.rm = T)
        end = max(times, na.rm = T)
      }
    }
    
    # figure out the times where we will do evaluations
    times<-seq(from = start, to=end,by = time.interval)
    
    tel<-as.data.frame.networkDynamic(nd)
    if (include.censored){  # detrmine of censored/truncated onset ties will be included in count
      formation<-sapply(times,function(t){sum(tel$onset==t)})
    } else {
      formation<-sapply(times,function(t){sum(tel$onset==t & !tel$onset.censored)})
    }
    
    if(result.type=='fraction'){
      # compute the number of empty dyads at each time point
      emptyDyads<-sapply(times,function(t){emptyDyadCount(nd,at=t)})
      # since we want to compare to number of empty dyads in 'previous' time point, add in the number of forming ties
      # (because ties that form now must have been empty in previous step)
      emptyDyads<-emptyDyads+formation
      formation<-formation/emptyDyads                 
    }
    return(ts(formation,start=start,end=times[length(times)],deltat=time.interval))
}

tEdgeDissolution<-function(nd, start, end, time.interval=1,result.type=c('count','fraction'),include.censored=FALSE){
  result.type=match.arg(result.type)
  if(missing(start) | missing(end)){
    times <- get.change.times(nd)
    if (length(times) == 0) {
      warning("network does not appear to have any time range information and start and end parameters not provided. Using start=0 end=0")
      start = 0
      end = 0
    } else {
      times[times == Inf] <- NA
      times[times == -Inf] <- NA
      start = min(times, na.rm = T)
      end = max(times, na.rm = T)
    }
  }
  
  # figure out the times where we will do evaluations
  times<-seq(from = start, to=end,by = time.interval)
  
  tel<-as.data.frame.networkDynamic(nd)
  if(include.censored){ # determine if terminus.censored ties will be included in the count
    dissolution<-sapply(times,function(t){sum(tel$terminus==t)})
  } else {
    dissolution<-sapply(times,function(t){sum(tel$terminus==t & !tel$terminus.censored)})
  }
  if(result.type=='fraction'){
    # compute the number of existing ties 
    activeECount<-sapply(times,function(t){network.edgecount.active(nd,at=t)})
    # add the number of disolving ties (since they would be counted as inactive above)
    # but recalculate to deal with spells where onset==terminus because they would be counted as active
    activeECount<-activeECount+sapply(times,function(t){sum(tel$terminus==t & tel$onset!=t)})
    dissolution<-dissolution/activeECount
  }
  return(ts(dissolution,start=start,end=times[length(times)],deltat=time.interval))
}

edgeFormationAt<-function(nd,at){
  tel<-as.data.frame.networkDynamic(nd)
  return(sum(tel$onset==at))
}

edgeDissolutionAt<-function(nd,at){
  tel<-as.data.frame.networkDynamic(nd)
  return(sum(tel$terminus==at))
}


# functions for calculating durations of verticse
vertexDuration<-function(nd,mode=c('duration','counts'),subject=c('vertices','spells'), v=seq_len(network.size(nd)), active.default=TRUE){
  # if network has no vertices, return nothing
  if (network.size(nd)==0){
    return(numeric(0))
  }
  
  del<-get.vertex.activity(nd,v=v,as.spellList = TRUE,active.default=active.default)
  # if looking per edge, group by edge id
  mode<-match.arg(mode)
  subject<-match.arg(subject)
  
  
  # determine function type for aggregation
  aggFun<-'sum' # function to use for aggregation sum = 'duration'
  if (mode=='counts'){
    aggFun<-'length'  # function to use to count events
  }
  
  # determine unit of aggragation
  if (subject=='vertices'){
    del<-aggregate.data.frame(del[,c('duration','vertex.id')],by=list(vertices=del$vertex.id),FUN=aggFun)
  }  else {
    # if looking at spells, just use the raw frame
    del<-aggregate.data.frame(del[,'duration',drop=FALSE],by=list(seq_len(nrow(del))),FUN=aggFun) 
  }
  
  return(del$duration)
  
}

# return the total amount of time that each vertex was tied via active edges
tiedDuration<-function(nd, mode=c('duration','counts'),active.default=TRUE,neighborhood=c('out','in','combined')){
  neighborhood<-match.arg(neighborhood)
  mode<-match.arg(mode)
  if(!is.directed(nd)){
    neighborhood<-'combined'
  }
  # determine function type for aggregation
  aggFun<-'sum' # function to use for aggregation sum = 'duration'
  if (mode=='counts'){
    aggFun<-'length'  # function to use to count events
  }
  bounds<-get_bounds(nd)
  spls<-as.data.frame.networkDynamic(nd,start=bounds[1],end=bounds[2],active.default=active.default)
  durations<-rep(0,network.size(nd))
  if(nrow(spls)>0){  # make sure there are some spells to aggregate
    if(neighborhood=='out'){
      connectDur<-aggregate(spls['duration'],list(spls$tail),aggFun)
      durations[connectDur[,1]]<-connectDur[,2]
    } else if (neighborhood=='in'){
      connectDur<-aggregate(spls['duration'],list(spls$head),aggFun)
      durations[connectDur[,1]]<-connectDur[,2]
    } else {  # ngh is combined
      connectDurHead<-aggregate(spls['duration'],list(spls$head),aggFun)
      connectDurTail<-aggregate(spls['duration'],list(spls$tail),aggFun)
      durations[connectDurHead[,1]]<-connectDurHead[,2]
      durations[connectDurTail[,1]]<-durations[connectDurTail[,1]]+connectDurTail[,2]
    }
  }
  return(durations)
}


# how much model clock time does it take on average for a single edge to change?
# divide the active vertex duration by ca ount the number of non-censored toggles in the network
meanTimeToChange<-function(nD){
  tel<-as.data.frame.networkDynamic(nD)
  bounds<-get_bounds(nD)
  # total number of changes is the number of on toggles (non-onset censored)
  #                               + number of off toggles (non-terminus censored)
  changeCount<-sum(!tel$onset.censored)+sum(!tel$terminus.censored)
  return((bounds[2]-bounds[1])*network.size(nD)/changeCount  )
}


nghOverlap<-function(nd, v, alter,neighborhood =c("out", "in", "combined"), type=c('bounds','total')){
  neighborhood<-match.arg(neighborhood)
  type<-match.arg(type)
  # is alter a neighbor of v?
  lapDur<-NA
  eids<-get.edgeIDs(nd,v=v,alter=alter,neighborhood = neighborhood)
  if (length(eids)>0){
    if (type=='bounds'){
      bounds<-range(get.edge.activity(nd,e = eids,as.spellList = TRUE)[c('onset','terminus')])
      lapDur<-bounds[2]-bounds[1]
    } else {
      # sum up all the spells on the edges
      lapDur<-sum(get.edge.activity(nd,e = eids,as.spellList = TRUE)$duration)
    }
  }
  return(lapDur)
}

# calculate the earliest time a path leaving v could reach alter (NOT earliest time leaving v), and the latest time a path leaving v could leave v (NOT latest time could arrive alter)
# I don't think this is a useful metric yet :-(
pathBounds<-function(nd, v, alter){
  netBounds<-get_bounds(nd)
  fwdDist<-paths.fwd.earliest(nd = nd,v=v,alter = alter,start=netBounds[1],end=netBounds[2])
  # check if there is any path at all
  if (is.infinite(fwdDist$tdist[alter])){
    return(c(NA,NA))
  }
  
  # find the neighbor on the first step of the early path from v leading to alter
  firstEarlyNgh<- NA
  prev<-alter
  while(prev != v){
    firstEarlyNgh<-fwdDist$previous[prev]
    prev<-fwdDist$previous[firstEarlyNgh]
  }
  bkwdDist<-paths.bkwd.latest(nd = nd,v=alter,alter=v,start=netBounds[1],end=netBounds[2])
  
  # find the neighbor on the last step of the late path from alter to v 
  # remember that path is backwards!
  lastLateNgh<-NA
  nxt<-v
  while(nxt != alter){
    lastLateNgh<-bkwdDist$previous[nxt]
    nxt<-bkwdDist$previous[lastLateNgh]
  }
  
  return(c(fwdDist$tdist[firstEarlyNgh]+fwdDist$start,bkwdDist$end-bkwdDist$tdist[lastLateNgh]))
}

# compute the number of empty dyads possible in a network
# allowing for network directedness, loops, 
emptyDyadCount <-function(net,at=NULL){
  if(is.multiplex(net)){
    stop("can not compute possible number of free dyads for multiplex networks")
  }
  n<-network.size(net)
  
  if(is.bipartite(net)){
    # need to calculate possible dyads based on feasible ties between the two modes
    firstSize<-net%n%'bipartite'
    secondSize<-n-firstSize
    if(has.loops(net)){
      # self loops don't make sense for bipartite networks
      warning('this network is bipartite and permits self-loop edges, loops ignored')
    }
    if(is.directed(net)){
      maxdyads<-firstSize*secondSize
    } else {
      maxdyads<-firstSize*secondSize/2
    }
  } else { # non- bipartite network
    if(is.directed(net)){
      maxdyads<-n*n
      # if it doesn't support loops, subtract the diagonal
      if(!has.loops(net)){
        maxdyads<- maxdyads-n
      }
    } else {
      maxdyads<- (n*(n-1))/2
      # if it has loops, add the diagnonal
      if(has.loops(net)){
        maxdyads<- maxdyads+n
      }
    }
  }
  
  # now subtract the number of existing edges
  if(is.networkDynamic(net)){
    if(is.null(at)){
      stop("networkDynamic object must include at parameter")
    }
    maxdyads<-maxdyads-network.edgecount.active(net,at=at,)
  } else {
    maxdyads<-maxdyads-network.edgecount(net)
  }
  
  return(maxdyads)
}
