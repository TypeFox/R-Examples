# draft version of a temporal mean degree measure
# assume input is in discrete time, 1 time step per unit
# TODO: include aggregate duration info
# TODO: how should degrees of non-active vertices be handled
tDegree<-function(nd, start, end, time.interval=1,cmode=c('freeman','indegree','outdegree')){
  
  if(missing(start) | missing(end)){
    times <- get.change.times(nd)
    if (length(times) == 0) {
      warning("network does not appear to have any dynamic information. Using start=0 end=1")
      start = 0
      end = 0
    }
    times[times == Inf] <- NA
    times[times == -Inf] <- NA
    start = min(times, na.rm = T)
    end = max(times, na.rm = T)
  }
  
  # figure out the times where we will do evaluations
  times<-seq(from = start, to=end,by = time.interval)
  
  cmode<-match.arg(cmode)
  # map cmode into a neighborhood type
  cmode<-switch(cmode,
         freeman='combined',
         indegree='in',
         outdegree='out')
  
  # allocate a matrix to hold the degrees of each vertex at each evaluation time
  degrees<-matrix(NA,ncol = network.size(nd),nrow=length(times))
  colnames(degrees) <- network.vertex.names(nd)
  
  # loop over times and compute degree of each vertex
  # TODO: would it be more efficient to stay on a single vertex and loop over time?
  for (i in seq_len(length(times))){
    at<-times[i]
    vids<-which(is.active(nd,at=at,v = seq_len(network.size(nd))))
    
    degreesAt<-sapply(vids,function(v){
      # get list of incident edges
      # TODO: need to correct for indegree, out degree ,etc
      eids<-get.edgeIDs(nd,v,neighborhood = cmode)
      # find out which ones are active
      active<-is.active(nd,e=eids,at = at)
      sum(active)
    })
    # deal with case of no active vertices where sapply returns list()
    if (!is.list(degreesAt)){
      degrees[i,vids]<-degreesAt
    }
    
  }
  # convert time timeseries object
  degrees<-ts(degrees,start=start,end=times[length(times)],deltat=time.interval)
  return(degrees)
}



