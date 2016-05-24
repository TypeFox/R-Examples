
######################################################################################################################################
######################################################################################################################################
### HiSSE Simulator -- work in progress 
######################################################################################################################################
######################################################################################################################################

#written by Brian O

SimulateHisse <- function(turnover.rates, eps.values, transition.rates, max.taxa=Inf, max.t=Inf, max.wall.time=Inf, x0, nstart=1, checkpoint.file=NULL, checkpoint.frequency=100, checkpoint.start.object=NULL, override.safeties=FALSE) {
    id <- living <- state <- NULL
    if(!is.finite(max.taxa) & !is.finite(max.t) & !is.finite(max.wall.time)) {
		if(!override.safeties) {
			stop("You have to limit the number of taxa, the tree height, and/or the actual run time. With current settings, hisse will grow a tree to infite size and height until the death of the universe. Or, until all the taxa in the simulation go extinct.")
		} else {
			warning("With current settings, hisse will grow a tree to infite size and height until the death of the universe. Or, until all the taxa in the simulation go extinct. Normally the program would throw an error, but you claim to know what you're doing (override.safeties==TRUE). We would strongly advise you to have checkpointing running. You may also want to buy a carbon offset for the CPU-years you might burn through during this simulation.")		
		}
	}
	start <- Sys.time()
	if(length(turnover.rates) != length(eps.values)) {
		stop("need to have same number of turnover and eps rates")	
	}
	if(length(turnover.rates) != dim(transition.rates)[1]) {
		stop("need to have same number of turnover and eps rates as rows in the transition matrix")	
	}
	if(length(turnover.rates) != dim(transition.rates)[1]) {
		stop("need to have same number of turnover and eps rates as columns in the transition matrix")	
	}
	state.levels <- c(0:(length(turnover.rates)-1))
	states <- factor(as.character(state.levels), levels=state.levels)
	results <- NA
	birth.rates <- GetBirthRate(turnover.rates, eps.values)
	death.rates <- GetDeathRate(turnover.rates, eps.values)
	diag(transition.rates) <- NA
	birth.counts <- 0*birth.rates
	death.counts <- 0*death.rates
	transition.counts <- 0*transition.rates

	if(!is.null(checkpoint.start.object)) {
		results <- checkpoint.start.object$results
		birth.counts <- checkpoint.start.object$birth.counts
		death.counts <- checkpoint.start.object$death.counts
		transition.counts <- checkpoint.start.object$transition.counts
	} else {
		results <- data.table(anc=NA, id=1, state=factor(x0, state.levels), length=0, height=0, living=TRUE, descendants=FALSE)
		for(additional.row in sequence(nstart-1)) {
			results <- rbind(results, data.table(anc=NA, id=additional.row+1, state=factor(x0, state.levels), length=0, height=0, living=TRUE, descendants=FALSE))
		}
	}
	setkey(results, id, living)
	keep.running <- TRUE
	rep.count <- 0
	while(keep.running) {
		rep.count <- rep.count+1
		tip.state.counts <- table(subset(results, living)$state)
		birth.rates.actual <- birth.rates * tip.state.counts
		death.rates.actual <- death.rates * tip.state.counts
		transition.rates.actual <- apply(transition.rates, 2, Multiply, y=tip.state.counts)
		birth.wait.times <- suppressWarnings(rexp(n=length(birth.rates.actual), birth.rates.actual))
		death.wait.times <- suppressWarnings(rexp(n=length(death.rates.actual), death.rates.actual))
		transition.wait.times <- suppressWarnings(matrix(rexp(n=length(transition.rates.actual), transition.rates.actual), nrow=dim(transition.rates.actual)[1])) #the NAs will be an issue
		min.times <- suppressWarnings(c(min(birth.wait.times, na.rm=TRUE), min(death.wait.times, na.rm=TRUE), min(transition.wait.times, na.rm=TRUE)))
		if((min(min.times, na.rm=TRUE)+max(subset(results, living)$height)) > max.t) { #gone too long
			keep.running <- FALSE	
			time.dif <- max.t - max(subset(results, living)$height)
			results[which(results$living),]$height <- results[which(results$living),]$height + time.dif
			results[which(results$living),]$length <- results[which(results$living),]$length + time.dif
		}
		if(keep.running) {
			results[which(results$living),]$height <- results[which(results$living),]$height + min(min.times)
			results[which(results$living),]$length <- results[which(results$living),]$length + min(min.times)

			if(which.min(min.times)==1) { #birth
				birth.counts[which.min(birth.wait.times)] <- birth.counts[which.min(birth.wait.times)]+1
				potential.lucky.taxa <- subset(results, living & state==states[which.min(birth.wait.times)])$id
				lucky.taxon <- potential.lucky.taxa[sample.int(length(potential.lucky.taxa), 1)]
				results[which(id==lucky.taxon),]$living <- FALSE
				results[which(id==lucky.taxon),]$descendants <- TRUE
				results <- rbind(results, data.table(anc=lucky.taxon, id=max(results$id)+1, state=subset(results, id==lucky.taxon)$state, length=0, height=subset(results, id==lucky.taxon)$height, living=TRUE, descendants=FALSE))
				results <- rbind(results, data.table(anc=lucky.taxon, id=max(results$id)+1, state=subset(results, id==lucky.taxon)$state, length=0, height=subset(results, id==lucky.taxon)$height, living=TRUE, descendants=FALSE))
			}
			if(which.min(min.times)==2) { #death
				death.counts[which.min(death.wait.times)] <- death.counts[which.min(death.wait.times)]+1
				potential.unlucky.taxa <- subset(results, living & state==states[which.min(death.wait.times)])$id
				unlucky.taxon <- potential.unlucky.taxa[sample.int(length(potential.unlucky.taxa), 1)]
				results[which(id==unlucky.taxon),]$living <- FALSE
			}
			if(which.min(min.times)==3) { #transition
				from.to <- which(transition.wait.times == min(transition.wait.times, na.rm=TRUE), arr.ind=TRUE)
				if (dim(from.to)[1] > 1) {
					from.to <- from.to[sample.int(dim(from.to)[1], 1),]	
				} else {
					from.to <- from.to[1,]	
				}
				transition.counts[from.to[1], from.to[2]] <- transition.counts[from.to[1], from.to[2]] +1
				transition.wait.times[from.to[1], from.to[2]] <- transition.wait.times[from.to[1], from.to[2]] + 1
				potential.changing.taxa <- subset(results, living & state==states[from.to[1]])$id
				changed.taxon <- potential.changing.taxa[sample.int(length(potential.changing.taxa), 1)]
				results[which(id==changed.taxon),]$state <- states[from.to[2]]
			}
			keep.running <- CheckKeepRunning(results, max.taxa, max.t, max.wall.time, start)
		}
		if(!is.null(checkpoint.file)) {
			if(rep.count %% checkpoint.frequency == 0) {
				checkpoint.result <- list(results=as.data.frame(results), birth.counts=birth.counts, death.counts=death.counts, transition.counts=transition.counts, n.surviving = dim(subset(results, living))[1])
				save(checkpoint.result, file=checkpoint.file)
			} 	
		}
	}
	return(list(results=as.data.frame(results), birth.counts=birth.counts, death.counts=death.counts, transition.counts=transition.counts, n.surviving = dim(subset(results, living))[1]))
}

Multiply <- function(x, y) { #I know, this is silly. It's like the joke about Wickham's addr package
	return(x*y)	
}

GetBirthRate <- function(turnover, eps) {
	return(turnover / (1+eps))	
}

GetDeathRate <- function(turnover, eps) {
	return((turnover * eps) / (1 + eps))	
}

CheckKeepRunning <- function(results, max.taxa=Inf, max.t=Inf, max.wall.time=Inf, start=NULL) {
    living <- NULL
    keep.running <- TRUE
	if(dim(subset(results, living))[1]<1) {
		keep.running <- FALSE	
	} 
	if(keep.running & dim(subset(results, living))[1]>=max.taxa) {
		keep.running <- FALSE
	}
	if(keep.running &  suppressWarnings(max(subset(results, living)$height)>=max.t)) {
		keep.running <- FALSE
	}
	if(keep.running & is.finite(max.wall.time)) { #could slow us down, so only check if needed
		if((as.numeric(Sys.time() - start)) > max.wall.time) {
			keep.running <- FALSE
		}
	}	
	return(keep.running)
}

SimToPhylo <- function(results, include.extinct=FALSE, drop.stem=TRUE) {
    living <- descendants <- NULL
    if(class(results)!="data.frame") {
		try(results <- 	results$results)
		if(class(results)!="data.frame") {
			stop("This requires a dataframe with simulation results")	
		}
	}
	if(dim(subset(results, living))[1]==0 & !include.extinct) {
		return(NA)
	}
	final.tips <- subset(results, !descendants)$id
	if(!include.extinct) {
		final.tips <- subset(results, living)$id
	}
	if(length(final.tips)==0) {
		return(NA)	
	}
	if(length(final.tips)==1) {
		return(structure(list(edge = structure(c(2L,1L), .Dim = c(1L, 2L)), edge.length = c(results[which(results$id==final.tips),]$height), tip.label = c("t1"), Nnode = 0L), .Names = c("edge", "edge.length", "tip.label", "Nnode"), class = "phylo"))
	}

	
	tips <- subset(results, !descendants)$id
	if( length(which(is.na(results$anc))) > 1) { #we do not have a stem, but start with node at base. Stick a stem on, then prune it off
		new.root.id <- min(results$anc, na.rm=TRUE)-1
		results[which(is.na(results$anc)),]$anc <- new.root.id
		results<- rbind(results[1,], results)
		results[1,]$anc <-NA
		results[1,]$id <- new.root.id
		results[1,]$length <- 0
		results[1,]$height <- 0		
		results[1,]$living <- FALSE
		results[1,]$descendants <- TRUE
		drop.stem <- TRUE
	}

	results$phylo.tipward.id <- NA
	results$phylo.tipward.id[which(!results$descendants)] <- sequence(length(tips))
	non.tips <- subset(results, descendants)$id
	results$phylo.tipward.id[which(results$descendants)] <- seq(from=(length(tips)+1), to=length(c(tips, non.tips)), by=1)
	results$phylo.rootward.id <- NA
	results$phylo.rootward.id <- sapply(results$anc, GetConversionOfAncestralNode, results=results)
	root.edge <- NULL
	if(drop.stem) {
		 results <- results[-which(is.na(results$anc)),]	
	} else {
		 root.edge <- results$length[which(is.na(results$anc))][1]
		 results <- results[-which(is.na(results$anc)),]	
	 }
	edge <-unname(cbind(as.numeric(results$phylo.rootward.id), as.numeric(results$phylo.tipward.id)))
	edge <- edge[order(edge[,2], decreasing=FALSE),]
	edge.length <- as.numeric(results$length[order(as.numeric(results$phylo.tipward.id), decreasing=FALSE)])
	Nnode <- length(non.tips)
	tip.label <- paste("t", sort(as.numeric(results$phylo.tipward.id[which(!results$descendants)]), decreasing=FALSE), sep="")
	tip.results <- results[which(!results$descendants),]
	states <- as.numeric(as.character(tip.results$state))
	names(states) <- paste("t", tip.results$phylo.tipward.id, sep="")
	phylo.return <- NA
	if(drop.stem) {
		phylo.return <- list(edge=edge, edge.length=edge.length, tip.label=tip.label, Nnode=Nnode)
	} else {
		phylo.return <- list(edge=edge, edge.length=edge.length, tip.label=tip.label, Nnode=Nnode, root.edge=root.edge)			
	}
	class(phylo.return)="phylo"
	phylo.return <- reorder(phylo.return)
	if(!include.extinct) {
		dead.tips <- subset(results, !descendants & !living)$phylo.tipward.id
		phylo.return <- drop.tip(phylo.return, dead.tips)	
		states <- states[!(names(states) %in% dead.tips)]
	}
	phylo.return$tip.state <- states
	return(phylo.return)
}

GetConversionOfAncestralNode <- function(anc.id, results) {
	return(results$phylo.tipward.id[which(results$id == anc.id)])
}