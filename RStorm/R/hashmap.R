# Functions to create and delete the local hashmap


SetHash <- function(name, data, ...){
	if(!is.list(RStorm.env$RStormData$hash)){
		RStorm.env$RStormData$hash <- list()
	}
  	RStorm.env$RStormData$hash[[name]] <- data
	TRUE
}

GetHash <- function(name, object=NULL){
	if(is.RStorm(object)){
		return(object$data[[name]])
	} else {
  		return(RStorm.env$RStormData$hash[[name]])
	}
}

.KillHash <- function(){
	RStorm.env$RStormData$hash <- NULL
}

GetHashNames <- function(object){
	if(!is.RStorm(object)){
		stop("Please provide an RStorm result object.")
	}
	print(names(object$data))
}

GetHashList <- function(object=NULL){
	if(is.RStorm(object)){
		return(object$data)
	} else {
  		return(RStorm.env$RStormData$hash)
	}
}
