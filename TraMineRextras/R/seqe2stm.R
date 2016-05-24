
seqe2stm <- function(events, dropMatrix=NULL, dropList=NULL, firstState="None"){
	nevent <- length(events)
	dropM <- NULL
	if(!is.null(dropList)){
		dropM <- matrix(FALSE, nrow=nevent,ncol=nevent, dimnames=list(events, events))
		for(i in 1:length(dropList)) {
			for(j in 1:length(dropList[[i]])) {
				## print(names(dropList)[i])
				## print(dropList[[i]][j])
				dropM[names(dropList)[i], dropList[[i]][j]] <- TRUE
			}
		}
	}
	else if(!is.null(dropMatrix)){
		if(all.equal(dim(dropMatrix),c(nevent, nevent) )){
			dropM <- dropMatrix
		}
		else {
			dropM <- matrix(FALSE, nrow=nevent, ncol=nevent)
		}
		dimnames(dropM) <- list(events, events)
		for(i in 1:nrow(dropMatrix)){
			dropM[rownames(dropMatrix)[i], colnames(dropMatrix)] <- dropMatrix[i,]
		}
	}
	nbstate <- 2^nevent
	stmatL <- matrix(TRUE, nrow=nbstate, ncol=nevent)
	colnames(stmatL) <- events
	stmat <- matrix("", nrow=nbstate, ncol=nevent)
	colnames(stmat) <- events
	stateNames <- function(occ, add=NULL){
		if(!is.null(add)){
			occ[add] <- TRUE
			if (!is.null(dropM)) {
				occ[dropM[add,]] <- FALSE
			}
		}
		if (sum(occ)==0) {
			return(firstState)
		}
		return(paste(events[occ], sep=".", collapse="."))
	}
	nr <- 1
	env <- environment()
	
	recursiveStateEvent <- function(i, occured){
		if(i>nevent) {
			env$stmatL[nr, ] <- occured
			env$nr <- nr+1
		} else {
			occured[i] <- FALSE
			recursiveStateEvent(i+1, occured)
			occured[i] <- TRUE
			recursiveStateEvent(i+1, occured)
		}
		invisible(TRUE)
	}
	recursiveStateEvent(1, rep(FALSE, nevent))
	stlist <- character(nbstate)
	for(i in 1:nbstate){
		stlist[i] <- stateNames(stmatL[i,])
		for(j in 1:nevent){
			occured <- stmatL[i,]
			stmat[i,j] <- stateNames(occured, add=j)
		}
	}
	rownames(stmat) <- stlist
	# rownamesMatrix <- matrix(rep(rownames(stmat), ncol(stmat)), nrow=nrow(stmat))
	# print(stmat)
	# print(rownamesMatrix)
	# print(stmat!=rownamesMatrix)
	# print(stmat[stmat!=rownamesMatrix])
	return(stmat[unique(c(firstState, stmat)),])
}