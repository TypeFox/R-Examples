get.cycle.matrix <-
function(dataset, CycleBreaks, rowLimit=NULL){
	limits <- NULL
  if(is.matrix(dataset)){
    if(is.null(rowLimit)){
  	for(i in 1:dim(CycleBreaks)[2]){
			breaks <- na.omit(CycleBreaks[,i])
			for(a in 1:(length(breaks)-1)){
				start <- breaks[a]
				end <- breaks[a+1]
				cycle <- dataset[start:end,i]
				limits <- append(limits, length(cycle))
				}
			}
		rowLimit <- max(limits)
		}
	cyclemat <- NULL	
	for(i in 1:dim(CycleBreaks)[2]){
		breaks <- na.omit(CycleBreaks[,i])	
		for(a in 1:(length(breaks)-1)){
			start <- breaks[a]
			end <- breaks[a+1]
			cycle <- dataset[start:end,i]
			cyclemat <- add.col(cyclemat, cycle, MaxTime=rowLimit) 
			}
		}
  }else{
  	if(is.null(rowLimit)){
  		breaks <- na.omit(CycleBreaks)
  		for(a in 1:(length(breaks)-1)){
  				start <- breaks[a]
  				end <- breaks[a+1]
  				cycle <- dataset[start:end]
  				limits <- append(limits, length(cycle))
  				}
  			
  		rowLimit <- max(limits)
  		}
  	cyclemat <- NULL	
  	
  		breaks <- na.omit(CycleBreaks)	
  		for(a in 1:(length(breaks)-1)){
  			start <- breaks[a]
  			end <- breaks[a+1]
  			cycle <- dataset[start:end]
  			cyclemat <- add.col(cyclemat, cycle, MaxTime=rowLimit) 
  			}
  		
  }
	return(cyclemat)
	}

