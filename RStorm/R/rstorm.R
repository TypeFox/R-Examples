# Dependencies
library(plyr)


# The RStrorm Environment:
RStorm.env <- new.env()
RStorm.env$RStormData <- list()

RStorm <- function(topology, 
					.verbose=TRUE, 
					.debug = FALSE, 
					.batches=100, ...){
	
	# Assign to the environment:
	RStorm.env$RStormData <- list("hash", "emit", "track")
	
	# Get stuff from topology
	spout <- topology$spout
	bolts <- topology$bolts
	finalize <- topology$finalize
	
	# Run through batches of the Spout:
	for(i in seq(1, nrow(spout), by=.batches)){
		
		#if(.verbose){print stuff}
		
		to <- ifelse((i+.batches)<nrow(spout), (i+.batches-1), nrow(spout))
		
		# Assign the first batch of rows to emit.0
		.BatchEmit(spout[i:to,,drop=FALSE], "emit.0");
		
		# Run each of the bolts
		for(j in c(1:length(bolts))){
			tmp.data <- .GetEmit(paste("emit.",bolts[[j]]$listen,sep=""))
	        alply(tmp.data, 1, bolts[[j]]$func, .name=paste("emit.",j,sep=""), boltID=bolts[[j]]$id, ...)
		}
		
		# Clean the emitted batches:
		if(.debug){  # If debug keep the last emitted set
			clean <- ifelse((i+.batches)<nrow(spout), TRUE, FALSE)
	    	if(clean){.RemoveEmit()}
		} else {  # else, remove everything
			.RemoveEmit()
		}
		
	}
	
	
	
	# Build the return object
	x <- list()
	x$data <- RStorm.env$RStormData$hash
	
	if(.debug){  # if in debug return emit sets
		cat("Returning debug emit sets...\n")
		x$debug <- RStorm.env$RStormData$emit
	}
	
	if(length(RStorm.env$RStormData$track) > 0){  # if data has been tracked add it
		x$track <- RStorm.env$RStormData$track
	}
	
	# Fun finalize function
	x$finalize <- finalize$func(RStorm.env$RStormData$hash)
	
	# assign class and return
	class(x) <- "RStorm"
	x
}


# Utilities


is.RStorm <- function(x){
  ifelse( class(x) == "RStorm", TRUE, FALSE )
}

print.RStorm <- function(x , ...){  ### EXTEND THIS ONE!!! ###
	cat("Name of the stored hashmaps:")
	print(names(x$data))
}

plot.RStorm <- function(x, ...){
	cat("There is currently no default plotting function available for RStorm objects.")
}




# Emit Function
Emit <- function(x, .name=NULL, ...){
	if(!is.Tuple(x)){
		stop("You can only emit Tuples. \nSee the examples in the documentation for more details.")
	}
	class(x) <- "data.frame"  # needed since RBind binds dataframes.
	if(is.data.frame(.GetEmit(.name))){
		x <- rbind(.GetEmit(.name), x)
	}
	RStorm.env$RStormData$emit[[.name]] <- data.frame(x)
	TRUE
}


# Internals
.GetEmit <- function(name){
	return(RStorm.env$RStormData$emit[[name]])
}

.BatchEmit <- function(x, name){
	RStorm.env$RStormData$emit[[name]] <- x
}

.RemoveEmit <- function(){
	RStorm.env$RStormData$emit <- NULL
}

