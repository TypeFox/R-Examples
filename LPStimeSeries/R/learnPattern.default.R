"learnPattern.default" <-
    function(x,
    	     segment.factor=c(0.05,0.95), 
	         random.seg=TRUE, target.diff=TRUE, segment.diff=TRUE, 
	         random.split=0,
             ntree=200,
             mtry=1,
             replace=FALSE,
             sampsize = if (replace) ceiling(0.632*nrow(x)) else nrow(x),
             maxdepth = 6,
             nodesize = 5,
	         do.trace=FALSE,
             keep.forest=TRUE,
             oob.pred=FALSE,
             keep.errors=FALSE, 
             keep.inbag=FALSE,   
             ...) {

	
    if(!is.matrix(x)){
		if(length(x)>0){ #single time series
			x <- t(as.matrix(x))
		}
		else{
			stop("data (x) has 0 rows")
		}   
    } 
     
    if (!is.numeric(x)) stop("data (x) is not numeric") 
        
    sampsize <- sampsize  
    n <- nrow(x)
    p <- ncol(x)
	
	if(length(segment.factor)>1){
		random.seg <- TRUE
		segment.factor <- sort(segment.factor)
	} else {
		random.seg <- FALSE
	}
    if (n == 0) stop("data (x) has 0 rows")
    
	ncat <- rep(1, p) #indicator for categorical variables (future purposes)
	
    ## overcome R's lazy evaluation:
    keep.forest <- keep.forest

    ## Make sure mtry is in reasonable range.
    if (mtry != 1)
        warning("invalid mtry: reset to within valid range")
    mtry <- 1
	
	if(oob.pred && target.diff){
		warning("Target segment cannot be difference series for OOB predictions, resetting target.diff=FALSE")
		target.diff=FALSE
	}
		
		
    ## Check for NAs.
    if (any(is.na(x))) stop("NA not permitted in predictors")

    ## Compiled code expects variables in rows and time series in columns.
    x <- t(x)
    storage.mode(x) <- "double"
    nt <- if (keep.forest) ntree else 1
    
    ## possible total # of nodes
    nrnodes <- 2^(maxdepth+1) - 1

	if(!replace){
		keep.inbag <- FALSE
		oob.pred <- FALSE
	} else {
		if(sampsize > n) stop("Sample size cannot be larger than number of time series")
		keep.inbag <- TRUE
		keep.errors <- TRUE
	}
		
	if(random.seg){
		segment.length <- runif(ntree,segment.factor[1],segment.factor[2])
	} else {
		segment.length <- rep(segment.factor[1],ntree)
	}

	rfout <- .C("regRF_time_series",
                    x,
                    as.double(segment.length),
                    as.integer(random.split),
                    as.integer(target.diff),
					as.integer(segment.diff),
                    as.integer(c(n, p)),
                    as.integer(sampsize),
                    as.integer(nodesize),
                    as.integer(nrnodes),
                    as.integer(ntree),
                    as.integer(mtry),
                    as.integer(ncat),
                    as.integer(do.trace),  
                    as.integer(oob.pred),  
                    as.integer(keep.errors),  
                    target = integer(ntree),           
                    target.type = integer(ntree),   
                    ndbigtree = integer(ntree),
                    nodedepth = matrix(integer(nrnodes * nt), ncol=nt),
                    nodestatus = matrix(integer(nrnodes * nt), ncol=nt),
                    splitType = matrix(integer(nrnodes * nt), ncol=nt),
                    leftDaughter = matrix(integer(nrnodes * nt), ncol=nt),
                    rightDaughter = matrix(integer(nrnodes * nt), ncol=nt),
                    nodepred = matrix(double(nrnodes * nt), ncol=nt),
                    bestvar = matrix(integer(nrnodes * nt), ncol=nt),
                    xbestsplit = matrix(double(nrnodes * nt), ncol=nt),
                    keep = as.integer(c(keep.forest, keep.inbag)),
                    replace = as.integer(replace),
                    oobpredictions = if (oob.pred)
                       double(n * p) else double(1),
                    ooberrors = if (oob.pred)
                       double(ntree) else double(1),
                    inbag = if (keep.inbag)
                       matrix(integer(n * ntree), n) else integer(1),
                    errors = if (keep.errors)
                       double(ntree) else double(1),
                    PACKAGE="LPStimeSeries")[c(16:32)]
        ## Format the forest component, if present.
        if (keep.forest) {
            max.nodes <- max(rfout$ndbigtree)
            rfout$nodestatus <-
                rfout$nodestatus[1:max.nodes, , drop=FALSE]
            rfout$nodedepth <-
                rfout$nodedepth[1:max.nodes, , drop=FALSE]                
            rfout$splitType <-
                rfout$splitType[1:max.nodes, , drop=FALSE]               
            rfout$bestvar <-
                rfout$bestvar[1:max.nodes, , drop=FALSE]
            rfout$nodepred <-
                rfout$nodepred[1:max.nodes, , drop=FALSE]
            rfout$xbestsplit <-
                rfout$xbestsplit[1:max.nodes, , drop=FALSE]
            rfout$leftDaughter <-
                rfout$leftDaughter[1:max.nodes, , drop=FALSE]
            rfout$rightDaughter <-
                rfout$rightDaughter[1:max.nodes, , drop=FALSE]
        }
        cl <- match.call()
        cl[[1]] <- as.name("learnPattern")

		if (oob.pred) rfout$oobpredictions[rfout$oobpredictions==-999]=NA
		
        out <- list(call = cl,
                    type = "regression",
                    random.split = random.split,
					segment.factor = segment.factor,
					segment.length = segment.length,
					nobs = floor(segment.length*p),
                    ntree = ntree,
                    maxdepth = maxdepth,
                    mtry = mtry,
                    target = rfout$target,
                    target.type = rfout$target.type,
                    forest = if (keep.forest)
                      c(rfout[c("ndbigtree", "nodedepth", "nodestatus", "splitType", "leftDaughter",
                              "rightDaughter", "nodepred", "bestvar",
                              "xbestsplit")],
                      list(nrnodes=max.nodes)) else NULL,
					oobpredictions= if (oob.pred)
                      t(matrix(rfout$oobpredictions, ncol=n)) else NULL,
					ooberrors= if (oob.pred)
                      rfout$ooberrors else NULL,                    
                    inbag = if (keep.inbag)
                      matrix(rfout$inbag, nrow(rfout$inbag),ntree) else NULL,
					errors = if (keep.errors) rfout$errors else NULL)
                           
    class(out) <- "learnPattern"
    return(out)
}
