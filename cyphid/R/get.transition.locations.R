get.transition.locations <-
function(cyclemat, close.cycle, type, CycleBreaks = NULL, closebreaks, frameLimit = NULL){
	if(is.null(frameLimit)){frameLimit <- 3}
	cycle <- NULL
	sequence <- NULL
	for(i in 1:dim(cyclemat)[2]){
		if(type=="FCSC"){
			dat <- cyclemat[1:close.cycle[i],i]
			}
		if(type=="SOFO"){
			cycleLength <- length(na.omit(cyclemat[,i]))
			dat <- cyclemat[close.cycle[i]:cycleLength,i]
			}
		if(is.dividing.valid(dat) == TRUE){
			out <- get.trans(dat)
			if(do.transitions.exist(out$transmin) == FALSE){
				cycle <- append(cycle, NA)
				}
			if(do.transitions.exist(out$transmin) == TRUE){
				nframes.btwn.close.and.transition <- abs(length(dat) - out$minframes[out$transmin])
				if(nframes.btwn.close.and.transition > frameLimit){
						if(type=="FCSC"){
							cycle <- append(cycle, out$minframes[out$transmin])
							}
						if(type=="SOFO"){
							cycle <- append(cycle, (out$minframes[out$transmin] + close.cycle[i])-1)
							}
						}else{
							cycle <- append(cycle, NA)
							}	
					}
			}
		if(is.dividing.valid(dat) == FALSE){
			cycle <- append(cycle, NA)
			}
		}
	
	if(!is.null(CycleBreaks)){
		breaks.vector <- get.breaks.vector(CycleBreaks)
		sequence <- cycle + breaks.vector
		cycleCounts <- get.cycle.counts(closebreaks)
		mat <- create.transition.location.matrix(sequence, cycleCounts)
		}
	return(list(cycle=cycle, sequence=sequence, seq.mat=mat))
	}

