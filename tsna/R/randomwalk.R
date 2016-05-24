# function to create a random forward walk on a network
# NOTE that this is not an ergodic process (once it heads down one branch, it might not be able to find its way back)

# TODO: include active.default
# TODO: option to force choosing discrete time units (round) or a minimum step size
paths.fwd.randomwalk<-function(nd,v,start,end,nsteps=1000){
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
  # compute time bounds
  changes<-NULL
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
    if(is.null(changes)){
      changes<-get.change.times(nd)
    }
    end<-max(changes)
  }
  
  # draw sample times that the walk will be evaluated at
  times<-sort(runif(n = nsteps,min=start,max=end))
  tdist<-rep(Inf,network.size(nd))
  previous<-rep(0,network.size(nd))
  tdist[v]<-0
  gsteps<-rep(Inf,network.size(nd))
  gsteps[v]<-0
  visitcounts<-rep(0,network.size(nd))
  visitcounts[v]<-1
  
  # walk the walk
  u<-v
  for (t in times){
    # find neighbors of u active at t
    ngs<-get.neighborhood.active(nd,v=u,type = 'out',at=t)
    if(length(ngs)>0){
      if(length(ngs)==1){
        w <-ngs
      } else {
        # select one at random
        w<-sample(ngs,1)
      }
      # update distances
      # since path may walk back on same vertex, only update if it hasn't been done before
      if (is.infinite(tdist[w])){
        tdist[w]<-t
        previous[w]<-u
        gsteps[w]<-gsteps[u]+1
      }
      visitcounts[w]<-visitcounts[w]+1
      u<-w
    } 
    
  }
  out<-list(tdist=tdist,
              previous=previous,
              gsteps=gsteps,
              start=start,
              end=end,
              direction='fwd',
              type='randomwalk',
              visitcounts=visitcounts)
  class(out)<-c('tPath',class(out))
  return(out)
}
