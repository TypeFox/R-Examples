get.cycle.breaks <-
function(dataset, window, rowLimit=NULL){
	# Get number of observations (i.e, frames) in each sequence
	if(is.matrix(dataset)){
    nobs <- dim(dataset)[1]
    cyclecount <- dim(dataset)[2]
	}else{
    nobs <- length(dataset)
    cyclecount <- 1
	}
  
	# Get max gape locations in the sequence
	CycleIndex <- get.peaks(-dataset, span = window)
	frames <- seq(1,nobs)
	if(is.null(rowLimit)){
		limits <- NULL
		for(i in 1:cyclecount){
      if(is.matrix(dataset)){
		    breaks <- frames[as.vector(na.omit(CycleIndex[,i]))]
      }else{
        breaks <- frames[CycleIndex]
      }
		  limits <- append(limits, length(breaks))
			}
		rowLimit <- max(limits)
		}
	CycleBreaks <- NULL
  NewIndex <- NULL
	for(i in 1:cyclecount){
		if(is.matrix(dataset)){
  	    #breaks <- frames[CycleIndex[,i]]
        IndexVector <- CycleIndex[, i]
        breaks <- frames[IndexVector]
        break.check <- check.for.sequence.endpoint.breaks(breaks, IndexVector)
        breaks <- break.check$breaks
        IndexVector <- break.check$IndexVector
        CycleBreaks <- add.col(CycleBreaks, breaks, rowLimit)
        NewIndex <- cbind(NewIndex, IndexVector)
        }else{
          CycleBreaks <- frames[CycleIndex]
        }
		#CycleBreaks <- add.col(CycleBreaks, breaks, rowLimit)
		}
	return(list(CycleBreaks=CycleBreaks, CycleIndex=NewIndex))
	}

