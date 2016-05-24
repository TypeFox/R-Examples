# Tracking Functions

TrackRow <- function(name, row){
	
	if(!is.list(RStorm.env$RStormData$track)){
		RStorm.env$RStormData$track <- list()
	}	
	
	if(nrow(row) != 1){
    	stop("Tracked rows should be a single row dataframe")
  	}
	
	if(is.data.frame(RStorm.env$RStormData$track[[name]])){
		RStorm.env$RStormData$track[[name]] <- rbind(RStorm.env$RStormData$track[[name]], row)
	} else {
		RStorm.env$RStormData$track[[name]] <- data.frame(row)
	}
        TRUE
}

GetTrackNames <- function(x){
	if(!is.RStorm(x)){
		stop("Please provide an RStorm result object.")
	}
	return(names(x$track))
}

GetTrack <- function(name, x){
	return(x$track[[name]])
}
