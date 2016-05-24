"ztransform" <- 
function(formula,data,family=gaussian) {
	if (missing(data)) {
		if(is(formula,"formula")) 
			data <- environment(formula)
		else  
			data <- environment()
#		wasdata <- 0
	} else {
		if (is(data,"gwaa.data")) {
			data <- data@phdata
		} 
		else if (!is(data,"data.frame")) {
			stop("data argument should be of gwaa.data or data.frame class")
		}
#		attach(data,pos=2,warn.conflicts=FALSE)
#		wasdata <- 1
	}
	
	if (is.character(family)) 
           family <- get(family, mode = "function", envir = parent.frame())
	if (is.function(family)) 
           family <- family()
	if (is.null(family$family)) {
           print(family)
           stop("'family' not recognized")
	}
	
	if ( is(try(formula,silent=TRUE),"try-error") ) { 
		formula <- data[[as(match.call()[["formula"]],"character")]] 
	}
	
	if (is(formula,"formula")) {
#		mf <- model.frame(formula,data,na.action=na.omit,drop.unused.levels=TRUE)
		mf <- model.frame(formula,data,na.action=na.pass,drop.unused.levels=TRUE)
		mids <- complete.cases(mf)
		mf <- mf[mids,]
		y <- model.response(mf)
		desmat <- model.matrix(formula,mf)
		lmf <- glm.fit(desmat,y,family=family)
#		if (wasdata) 
#			mids <- rownames(data) %in% rownames(mf)
#		else 
		resid <- lmf$resid
#		print(formula)
	} else if (is(formula,"numeric") || is(formula,"integer") || is(formula,"double")) {
		y <- formula
		mids <- (!is.na(y))
		y <- y[mids]
		resid <- y
		if (length(unique(resid))==1) stop("trait is monomorphic")
		if (length(unique(resid))==2) stop("trait is binary")
	} else {
		stop("formula argument must be a formula or one of (numeric, integer, double)")
	}
	y <- (resid-mean(resid))/sd(resid)
#	if (wasdata==1) detach(data)
	tmeas <- as.logical(mids)
	out <- rep(NA,length(mids))
	out[tmeas] <- y
	out
}

