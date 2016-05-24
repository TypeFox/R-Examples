# Modified: 6 Sept 2015

collapseClips <-
function(
  rec,          # The wave object (or file) to be collapsed
  start.times,  # The start times of clips (s)
  end.times,    # The end times of clips (s)
  return.times=FALSE  # Set to TRUE to return the times
) {

  # If rec isn't already a wave object, but it is character, it is assumed to be a file path and file is read in. 
  rec <- getClip(rec, output="Wave")
  start.times[start.times<0] <- 0
  end.times[end.times>length(rec@left)/rec@samp.rate] <- length(rec@left)/rec@samp.rate
  times <- cbind(start.times, end.times)
  duration <- end.times - start.times
  cum.times <- data.frame(start.time.collapse=cumsum(duration)-duration, end.time.collapse=cumsum(duration))

  # Clip out pieces of rec, and save as a list
  waves <- apply(times, 1,function(x) cutWave(rec, from=x[1], to=x[2]))

  # Then combine them all with bind, using do.call so waves object can be used as the argument list for the ... argument in bind
  collapsed <- do.call(tuneR::bind, waves)
  if(return.times) return(list(wave=collapsed, times=cum.times)) else return(collapsed)
}
