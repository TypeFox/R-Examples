######################################################################################################################################
######################################################################################################################################
### TransMatMaker -- Builds transition rate matrix for easy use in the main function
######################################################################################################################################
######################################################################################################################################

TransMatMaker <- function(hidden.states=FALSE){
	if(hidden.states == FALSE){
		rate.mat <- matrix(NA, 2, 2)
		diag(rate.mat) <- 3
		rate.mat[is.na(rate.mat)] <- 1:2
		diag(rate.mat) <- NA
		rownames(rate.mat) <- c("(0)","(1)")
		colnames(rate.mat) <- c("(0)","(1)")					
	}else{
		rate.mat <- matrix(NA, 4, 4)
		diag(rate.mat) <- 13
		rate.mat[is.na(rate.mat)] <- 1:12
		diag(rate.mat) <- NA
		rownames(rate.mat) <- c("(0A)","(1A)","(0B)","(1B)")
		colnames(rate.mat) <- c("(0A)","(1A)","(0B)","(1B)")			
	}
	return(rate.mat)
}


######################################################################################################################################
######################################################################################################################################
### Various functions for dropping and setting equal parameters in a transition matrix.
######################################################################################################################################
######################################################################################################################################

ParDrop <- function(rate.mat.index=NULL, drop.par=NULL){
	if(is.null(rate.mat.index)){
		cat("Rate matrix needed.  See TransMatMaker() to create one.\n")
		return
	}
	if(is.null(drop.par)){
		cat("No parameters indicated to drop.  Original matrix returned.\n")
		return(rate.mat.index)
	}
	if(max(rate.mat.index,na.rm=TRUE) < max(drop.par,na.rm=TRUE)){
		cat("Some parameters selected for dropping were not in the original matrix.\n")
	}
	drop.par <- unique(drop.par) # in case parameters listed more than once in drop vector
	drop.par <- drop.par[order(drop.par)]
	max <- max(rate.mat.index,na.rm=TRUE)
	for(drop.which in 1:length(drop.par)){
		drop.locs <- which(rate.mat.index == drop.par[drop.which],arr.ind=TRUE)
		rate.mat.index[drop.locs] <- NA
	}
	max <- max - length(drop.par)
	exclude <- which(is.na(rate.mat.index))
	rate.mat.index[-exclude] <- 1:max
	rate.mat.index[is.na(rate.mat.index)] = 0
	diag(rate.mat.index) = NA
	return(rate.mat.index)
}


ParEqual <- function(rate.mat.index=NULL, eq.par=NULL){
	if(is.null(rate.mat.index)){
		cat("Rate matrix needed.  See TransMatMaker() to create one.\n")
		return
	}
	if(is.null(drop) || length(eq.par) < 2){
		cat("Fewer than two parameters indicated to equalize. Original matrix returned.\n")
		return(rate.mat.index)
	}
	too.big <- which(eq.par > max(rate.mat.index,na.rm=TRUE))
	if(length(too.big) > 0){
		cat("Some parameters selected for equalizing were not in the original matrix:\n")
		cat("Not in original rate.mat.index:",eq.par[too.big],"\n")
		cat("Original matrix returned.\n")
		return(rate.mat.index)
	}
	eq.par <- unique(eq.par)
	eq.par <- eq.par[order(eq.par)]	
	min <- min(eq.par) # rm.na unnecessary?
	
	#the decrement index will hold counters to decrement rate index
	dec.index <- matrix(0,length(rate.mat.index[,1]),length(rate.mat.index[1,]))
	for(eq.which in 2:length(eq.par)){
		to.eq <- which(rate.mat.index == eq.par[eq.which],arr.ind=TRUE)
		rate.mat.index[to.eq] <- min
	}
	#the decrement index will hold counters to decrement rate index
	dec.index <- matrix(0,length(rate.mat.index[,1]),length(rate.mat.index[1,]))
	for(eq.which in 2:length(eq.par)){
		to.dec <- which(rate.mat.index > eq.par[eq.which],arr.ind=TRUE) #greater than current decrementer
		dec.index[to.dec] <- dec.index[to.dec] + 1
	}
	rate.mat.index <- rate.mat.index - dec.index
	rate.mat.index[is.na(rate.mat.index)] = 0
	diag(rate.mat.index) = NA
	return(rate.mat.index)
}

