#  Part of the statnet package, http://statnetproject.org
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) in
#    http://statnetproject.org/attribution
#
#  Copyright 2014 the statnet development team
######################################################################

# functions for computing temporal density

# if regular graph density is the fraction of possible ties that exist
# temporal density would be the fraction of total possible tie-time that is active
# OR  fraction of time existing ties are active?


# general wrapper function
tEdgeDensity<-function(nd,mode=c('duration','event'),agg.unit=c('edge','dyad'),active.default=TRUE){
  if (!is.networkDynamic(nd)){
    stop('tEdgeDensity requires that its first argument is a networkDynamic object')
  }
  mode<-match.arg(mode)
  agg.unit<-match.arg(agg.unit)
  bounds<-get_bounds(nd)
  # trap some edge cases
  if (network.edgecount(nd)==0 | network.size(nd)==0){
    return(0)
  }
  
  # if we are doing durations..
  if (mode=='duration'){
    spls<-as.data.frame.networkDynamic(nd,start=bounds[1],end=bounds[2],active.default=active.default)
    total_dur<-sum(spls$duration)
    if (agg.unit=='edge'){ # at the edge level
      return(total_dur/(network.edgecount(nd)*(bounds[2]-bounds[1])))
    }
    if (agg.unit=='dyad'){ # at the dyad level
      # total number of possible dyads
      num_dyads<-network.dyadcount(nd)
      return(total_dur/(num_dyads*(bounds[2]-bounds[1])))
    }
    
  } else if (mode=='event') {  # we are doing event counts
      if (agg.unit=='edge'){
      spls<-get.edge.activity(nd)
      nulls<-sapply(spls, is.null)
      counts<-sapply(spls[!nulls], nrow)
      numSpells<-sum(counts)
      # TODO: don't count always active or always inactive because they don't change
      return(numSpells/(network.edgecount(nd)*bounds[2]-bounds[1]))
    } 
  }
  # if we get here, something is wrong
  stop('tEdgeDensity is not yet implemented for ',agg.unit,' ',mode,'.')
}


# TODO: should self-loops be counted? Only if loops=TRUE?
# TODO: what about multiplex networks?




# how many events are there per time step within the graph time period
edge_event_density<-function(nd){
  
  # trap some edge cases
  if (network.edgecount(nd)==0 | network.size(nd)==0){
    return(0)
  }
  bounds<-get_bounds(nd)
  # count the number of events occuring
  spls<-get.edge.activity(nd)
  nulls<-sapply(spls, is.null)
  counts<-sapply(spls[!nulls], nrow)
  numSpells<-sum(counts)
  # TODO: don't count always active or always inactive because they don't change
  return(numSpells/(network.edgecount(nd)*bounds[2]-bounds[1]))
}



dyad_duration_density<-function(nd,active.default=TRUE){
  # trap some edge cases
  if (network.edgecount(nd)==0 | network.size(nd)==0){
    return(0)
  }
  bounds<-get_bounds(nd)
  spls<-as.data.frame.networkDynamic(nd,start=bounds[1],end=bounds[2],active.default=active.default)
  total_dur<-sum(spls$duration)
  # total number of possible dyads
  num_dyads<-network.dyadcount(nd)  # corrects for self loops as of network 1.13

  return(total_dur/(num_dyads*(bounds[2]-bounds[1])))
}

edge_duration_density<-function(nd,active.default=TRUE){
  # trap some edge cases
  if (network.edgecount(nd)==0 | network.size(nd)==0){
    return(NA)
  }
  bounds<-get_bounds(nd)
  spls<-as.data.frame.networkDynamic(nd,start=bounds[1],end=bounds[2],active.default=active.default)
  total_dur<-sum(spls$duration)
  return(total_dur/(network.edgecount(nd)*(bounds[2]-bounds[1])))
}


