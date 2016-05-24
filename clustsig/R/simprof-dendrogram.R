simprof.plot <- function(results, leafcolors=NA, plot=TRUE, fill=TRUE, leaflab="perpendicular", siglinetype=1){
	firstfind <- rep(TRUE, results$numgroups)
	colour <- local({
		colorizer <- function(i, siggroups, leafcolors, fill, siglinetype){
			#
			# i = the current node (passed by dendrapply I assume)
			#
			# siggroups should be $significantclusters from the results of simprof
			# so siggroups is a list of vectors
			#
			# leafcolors should be a vector of colors equal in number to the number length of siggroups
			# leafcolors <- rainbow(length(siggroups)) basically
			#
			
			props <- attributes(i)
			if(is.leaf(i)){
				leafID <- props$label
				group <- findGroup(leafID, siggroups)
				color <- leafcolors[group]
				attr(i, "edgePar") <- c(props$edgePar, list(col=c(color, color), lty=siglinetype[group]))
				}
			if(!is.leaf(i) && fill){
				mems <- labels(i)
				unif <- pure(mems, siggroups)
				if (unif > 0){
					if (firstfind[unif]){
						firstfind[unif] <<- FALSE
						}
					else {
						color <- leafcolors[unif]
						attr(i, "edgePar") <- c(props$edgePar, list(col=c(color,color), lty=siglinetype[unif]))
						}
					}
				}
			i
			}
		})
	siggroups <- results$significantclusters
	if (typeof(siglinetype) == "character"){ # not && to get rid of warning
		if (siglinetype == "different"){
			siglinetype <- rev(c(1:6)) # make it so that "solid" is used last.
			if (length(siggroups) > length(siglinetype)){
				warning("More significant groups than there are valid line types. Some line types will be repeated.")
				while (length(siggroups) > length(siglinetype))
					siglinetype <- append(siglinetype, siglinetype) # just keep adding on until we are over how many we need
				}
			}
		else 
			siglinetype <- rep(1, length(siggroups)) # something went wrong
		}
	else if (length(siglinetype) == 1 && siglinetype[1] >= 1 && siglinetype[1] <= 6)
		siglinetype <- rep(siglinetype, length(siggroups))
	else if (length(siggroups) > length(siglinetype)){
		warning("Fewer line types supplied than there are significant groups.")
		if (length(siglinetype) < 6){
			siglinetype <- rev(c(1:6))
			if (length(siggroups) > length(siglinetype))
				while (length(siggroups) > length(siglinetype))
					siglinetype <- append(siglinetype, siglinetype) # just keep adding on until we are over how many we need
			}
		}
	
	if (is.na(leafcolors[1]))
		leafcolors <- rainbow(length(results$significantclusters))
	else if (length(leafcolors) == 1)
		leafcolors <- rep(leafcolors[1], length(siggroups))
		
	
	dend <- dendrapply(as.dendrogram(results$hclust), colour, siggroups=results$significantclusters, leafcolors=leafcolors, fill=fill, siglinetype=siglinetype)
	if(plot)
		plot(dend, leaflab=leaflab)
	return(dend)
	}
	
findGroup <- function(leafID, siggroups){
	for (j in 1:length(siggroups)){
		if (length(which(siggroups[[j]] == leafID) != 0)){
			return(j) # this should be sufficient
			}
		}
	}

pure <- function(mems, siggroups){
	for (j in 1:length(siggroups))
		if (all(!is.na(match(mems, siggroups[[j]]))))
			return(j) # return the group they all belong to
	return(-1)
	}
