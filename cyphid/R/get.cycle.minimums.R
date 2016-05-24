get.cycle.minimums <-
function(dataset, CycleBreaks, cyclemat, rowLimit=NULL){
	# Get cycle minimums with respect to the cycle
	close.cycle <- NULL
	for(i in 1:dim(cyclemat)[2]){
		close.cycle <- append(close.cycle, which.max(cyclemat[,i]))
		}
	# Get cycle minimums with respect to the cycle sequence
	if(is.null(rowLimit)){
		limits <- NULL
		for(i in 1:dim(CycleBreaks)[2]){
			set1 <- na.omit(CycleBreaks[,i])
			ncycles <- length(set1)-1
			limits <- append(limits, ncycles)
			}
		rowLimit <- max(limits)
		}
	closebreaks <- NULL
	start <- 1
	for(i in 1:dim(CycleBreaks)[2]){
		set1 <- na.omit(CycleBreaks[,i])
		nbreaks <- length(set1)
		ncycles <- nbreaks-1
		open <- CycleBreaks[1:ncycles, i]
		end <- start + ncycles - 1
		close <- close.cycle[start:end]
		set2 <- open + close
		closebreaks <- add.col(closebreaks, set2, rowLimit)
		start <- end + 1
		}
	return(list(closebreaks=closebreaks, close.cycle=close.cycle))
	}

