## Create and add to a topology

Topology <- function(spout, name=NULL, .verbose=TRUE){
	x <- list()
	if(is.data.frame(spout)){
		x$spout <- spout
		lines <- dim(spout)[1]
		if(.verbose){
			cat(paste("Created a topology with a spout containing ", lines, "rows."))
		}
	} else {
		stop("Please provide a dataframe as a spout")
	}
	x$bolts <- list()
	x$finalize <- Bolt(function(x){FALSE})
	class(x) <- "Topology"
	x
}

is.Topology <- function(x){
  ifelse( class(x) == "Topology", TRUE, FALSE )
}

AddBolt <- function(topology, bolt, .verbose=TRUE){
	if(!is.Topology(topology)){
		stop("Please provide a topology")
	}
	if(!is.Bolt(bolt)){
		stop("Please provide a bolt")
	}
	place <- length(topology$bolts) + 1
 	topology$bolts[[place]] <- bolt
	if(.verbose){
		print(paste("Added bolt",bolt$name,"to position",place,"which listens to",bolt$listen))
	}
	topology
}



AddFinalize <- function(topology, bolt, .verbose=TRUE){
	topology$finalize <-  bolt
	cat("Added finalize function",bolt$name);
	return(topology)
}

ShowFinalize <- function(x){
	print(x$finalize)
}

ChangeSpout <- function(topology, spout){
	topology$spout <- spout
	topology
}


print.Topology <- function(x, ...){
	lines <- dim(x$spout)[1]
	cat(paste("Topology with a spout containing", lines, "rows","\n"))
	if(length(x$bolts) > 0){
		for(i in 1:length(x$bolts)){
			cat(paste(" - Bolt (",i,"): *", x$bolts[[i]]$name, "* listens to",  x$bolts[[i]]$listen,"\n"))
		}
	} else {
		cat("No bolts specified \n")
	}
	
	if((x$finalize$name[[1]]) == "function"){
	       cat("No finalize function specified \n")
	} else {
              cat("With finalize function",x$finalize$name,"\n")
	}

}

plot.Topology <- function(x, ...){
	cat("There is currently no default plotting function available for Topology objects.")
}


