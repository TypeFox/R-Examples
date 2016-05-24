# imat2 <- function(f, all = TRUE) {
  # f <- as.factor(f)
  # lev <- levels(f)
  # if(!all) lev <- lev[-length(lev)]
  # sapply(lev, "==", f) + 0
# }

# idat2 <- function(x, all = FALSE, keep = if("Freq" %in% names(x)) "Freq" else NULL) {
  # force(all)
  # X <- do.call(cbind, 
         # lapply(x, function(f) 
           # if(is.factor(f)) imat2(f, all) else NULL))
  # cbind(data.frame(X), Freq = x[, keep])
# }



imat <- function(x, allcat = TRUE) {
	#if(!is.null(dim(x))){
	#	if(!length(dim(x)==1)==1){
	#		stop("For dataframes and matrices please use idat instead.")
	#	}
	#	dim(x) <- NULL
	#}
  x <- as.factor(x)
  
  	if(any(is.na(x))){
  		levels(x) <- c(levels(x),"N/A")
		x[is.na(x)] <- "N/A"
	}
  
  lev <- levels(x)
  m <- length(lev)
  
  if(m == 1){
  	return(x)
  }

 	X <- matrix(0, n <- length(x), m <- length(lev) )
	X[(1:n) + n*(unclass(x)-1)] <- 1
	dimnames(X) <- list(names(x), lev)
		
	if( !allcat ){
  		return(X[,-m,drop=FALSE])
  	}
	
	return(X)
}

idat <- function(x, allcat = FALSE, keep = "Freq") {
 
 # convert varnames to indices
  if(is.character(keep[1])){
  	keep <- names(x) %in% keep
  	keep <- if(any(keep)) which(keep) else NULL
  }
  
  #if(! all(is.integer(keep)) ){
  #	stop("Argument keep must contain either indices or variable names.")
  #}
  if(is.null(keep)){
  	 X <- lapply(x, function(f) 
           	imat(f, allcat))
     labs <- names(x)
  }else{
  	 X <- lapply(x[,-keep], function(f) 
           	imat(f, allcat))
     labs <- names(x)[-keep]
  }
 
  nlvl <- sapply(X, ncol)
  X <- do.call(cbind, X)
  
  ret <- data.frame( cbind(X, x[, keep]))	
    names(ret) <- c(paste(  rep(labs ,nlvl), dimnames(X)[[2]],sep = ":"), names(x)[keep])
   
  return(ret)
}