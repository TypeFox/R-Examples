simprepTemporal <-
function(totaltime, regions=NULL, onsets, durations, TR, effectsize, accuracy=0.1, hrf=c("gamma","double-gamma","Balloon"), param=NULL){

	if(missing(hrf)){
		hrf <- "double-gamma"
	}
	if(!is.list(onsets)){
	  if(is.numeric(onsets)){
	    onsets <- list(onsets)
	  }else{
	    stop("Wrong format of onsets. See ?simprepTemporal for further details.")
	  }
	}
	if(!is.list(durations)){
	  if(is.numeric(durations)){
	    durations <- list(durations)
	  }else{
	    stop("Wrong format of durations. See ?simprepTemporal for further details.")
	  }
	}
	if(!is.list(effectsize)){
	  if(is.numeric(effectsize)){
	    effectsize <- list(effectsize)
	  }else{
	    stop("Wrong format of durations. See ?simprepTemporal for further details.")
	  }
	}
	if(!is.null(regions)){
	    if(regions!=1){
		if((length(onsets)!=regions)||(length(durations)!=regions)){
		  stop("Mismatch between number of regions and onsets and/or durations. See ?simprepTemporal for further details.")
		}
	    }
	} else {
		regions <- 1
	}
	if(!is.null(param)){
	    if(regions!=1){
	      if(length(param)!=regions){
		stop("Mismatch between number of regions and param. See ?simprepTemporal for further details.")
	      }
	    }
	}
	
	temp <- list()
	if(regions == 1){
		temp[[1]] <- list()
		name <- "region1"
		temp[[1]]$name <- name
		temp[[1]]$onsets <- onsets
		temp[[1]]$durations <- durations
                temp[[1]]$totaltime <- totaltime
		temp[[1]]$effectsize <- effectsize
		temp[[1]]$TR <- TR
		temp[[1]]$accuracy <- accuracy
		temp[[1]]$hrf <- hrf
		if(!is.null(param)){
		  temp[[1]]$param <- param
		} else {
		  temp[[1]]$param <- NULL
		}
	} else {
		for(i in 1:regions){
			temp[[i]] <- list()
			name <- paste("region", i, sep="")
			temp[[i]]$name <- name
			temp[[i]]$onsets <- onsets[[i]]
			temp[[i]]$durations <- durations[[i]]
                        temp[[i]]$totaltime <-totaltime
			temp[[i]]$effectsize <- effectsize[[i]]
			temp[[i]]$TR <-TR
			temp[[i]]$accuracy  <- accuracy
			temp[[i]]$hrf <- hrf
			if(!is.null(param)){
				temp[[i]]$param <- param[[i]]	
			} else {
				temp[[i]]$param <- NULL
			}
		}
	}
	return(temp)
}

