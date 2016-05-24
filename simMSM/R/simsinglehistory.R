simsinglehistory <- function(first.entry = 0, first.from, max.time, change.times, mpl, x.i, 
                             partial.markov.x, partial.markov.eta){
  current.exit <- current.entry <- first.entry
  current.to <- current.from <- first.from
  history.i <- pmX <- pmX.i <- NULL
  p <- length(x.i)
  while((current.entry < max.time) & (!is.null(mpl[[current.from]]$all.to))){
    eta.ij <- vector("list", length(mpl[[current.from]]$eta))
    for(j in 1:length(eta.ij)){
      eta.ij[[j]] <- mpl[[current.from]]$eta[[j]]
    }
    ## partial markov models: generate x!
    if(!is.null(partial.markov.x)){
      if(is.null(history.i)){
        ## pmX.i <- partial.markov.x(rbind(NULL, rep(0, 5)))
        pmX.i <- partial.markov.x(0, 0, 0, 0, 0)
        pmX <- pmX.i <- pmX.i * 0
      }else{
        pmX.i <- partial.markov.x(entry = history.i[, 1], exit = history.i[, 2], 
                                  from = history.i[, 3], to = history.i[, 4], 
                                  delta = history.i[, 5])
        pmX <- rbind(pmX, pmX.i)
      }
    }
    ## partial markov models: calculate partial markov etas!
    pme <- rep(0, length(mpl))
    if(!is.null(pmX.i)){
      for(j in 1:length(pme)){
        helpfun <- partial.markov.eta[[current.from]][[j]]
        pme[j] <- helpfun(pmX.i)
      }
    }
    ## sim exit and to:
    current.sim <- simto(current.entry, current.from, mpl, eta.ij = eta.ij, 
                         x.i = x.i, max.time = max.time, pme = exp(pme))
    delta.ij <- 1
    if(current.sim$exit.ij == current.entry){
      current.sim$exit.ij <- current.sim$exit.ij + 1e-10
    }
    current.entry <- current.sim$exit.ij
    current.from <- current.sim$to.ij
    if(!is.null(change.times)){
      for(cti in 1:length(change.times)){
        if((current.sim$entry.ij < change.times[cti]) & (current.sim$exit.ij > change.times[cti])){
          current.sim$exit.ij <- change.times[cti]
          delta.ij <- 0
          current.entry <- change.times[cti]
          current.from <- current.sim$from.ij
        }
      }
    }
    history.ij <- c(current.sim$entry.ij, current.sim$exit.ij, 
                    current.sim$from.ij, current.sim$to.ij, delta.ij)
    ## if observation is censored (due to change.time!):
    #if(history.ij[5] == 0){
    #  all.to.i <- mpl[[history.ij[3]]]$all.to
    #  if(length(all.to.i) > 1.5){
    #    for(i.plus in 2:length(all.to.i)){
    #      history.ij <- rbind(history.ij, history.ij)
    #      pmX <- rbind(pmX, pmX.i)
    #    }
    #    ## overwrite to state with mpl[[from.i]]$all.to:
    #    history.ij[, 4] <- all.to.i
    #  }
    #}
    history.i <- rbind(history.i, history.ij)
  }
  return(list(history.i = history.i, pmX = pmX))}