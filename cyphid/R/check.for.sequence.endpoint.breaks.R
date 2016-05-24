check.for.sequence.endpoint.breaks <-
function(breaks, IndexVector){
  FrameCount <- length(na.omit(IndexVector))
  breaks <- as.vector(na.omit(breaks))
  if(breaks[1]==1){
    end <- length(breaks)
    breaks <- breaks[2:end]
    IndexVector[1] <- FALSE
  }
  if(breaks[length(breaks)]==FrameCount){
    end <- length(breaks) - 1
    breaks <- breaks[1:end]
    IndexVector[FrameCount] <- FALSE
  }
  return(list(breaks=breaks, IndexVector=IndexVector)) 
}

