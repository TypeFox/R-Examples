# function to compute the rate at which new vertices can be reached for a network
#get a list of seed vertices
#for each seed
#walk the forward reachable path
#record the size of the forward reachable set
#divide by the duration of the network
reachableRate<-function(nD, start, end, seeds){
  if (missing(start) | missing(end)) {
    times <- get.change.times(nD)
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
  
  reachableRates<-sapply(seeds,function(s){
    reachableN<-sum(tPath(nD,v=s,start=start,end=end)$tdist<Inf)
    reachableN/(end-start)
  })
  #average over all the seeds
  return(mean(reachableRates))
}


# martina points out that this measure above may suffer from the boundry condition created by network size: As time goes on, there are fewer and fewer unreached vertices for the path to discover, so its rate will slow.  Perhaps a solution is to measure the paths from sets of random start times? 




# for a set of seeds, compute forward reachable path and times when reach occurs
# return a vector of 10 times in which each element is the mean number of vertices reached 
# (from seeds) at that time point
meanReachTimes<-function(nD, start, end, seeds){
  if (missing(start) | missing(end)) {
    times <- get.change.times(nD)
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
  distances<-lapply(seeds, function(s){
    tPath(nD,v =s,start = start,end=end,graph.step.time = 1 )$tdist
    })
  times<-seq(from=start,to=end,length.out = 10)
  means<-sapply(times,function(t){
    # compute the
    dAtT<-sapply(distances,function(d){
      sum(d<=t)
    })
    mean(dAtT)
  })
  names(means)<-times
  return(means)
  
}

# for a set of seets, compute forward reachable path and times when reach occurs
# return the times at which each seed reached at least num.targets vertices.

# TODO: should be able to do this faster with a customized tPath function 
# that stops when the required number of vertices are found
timeToReach<-function(nD, num.targets=round(network.size(nD)/2), start, end, seeds){
  if (missing(start) | missing(end)) {
    times <- get.change.times(nD)
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
  distances<-lapply(seeds, function(s){
    tPath(nD,v =s,start = start,end=end,graph.step.time = 1 )$tdist
  })

  reachTimes<-sapply(distances,function(d){
    d<-sort(d) # put reachable times in temporal order
    # find the time value at num.targets
    # (will be Inf if never reached that many)
    d[num.targets]
  })
  return(reachTimes)
}


