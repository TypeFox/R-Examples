mnl.spec <- function(f, data, choices=NULL, base.choice=1, varying = NULL, sep='.') {
	# Create an MNL specification object 
	
	f.char <- deparse(f)
	if (length(f.char) > 1) f.char <- paste(f.char, collapse='')
	f.spl <- unlist(strsplit(f.char, '|', fixed=TRUE))
	# everything before '|' means the same coefficient for all alternatives
	f1 <- formula(f.spl[1]) 
	tms <- terms(formula(f1), data = data)
    tms.labels <- attr(tms, "term.labels")
    tms.intercept <- attr(tms, "intercept")
    		
    tms.v.labels <- NULL
    if(length(f.spl) > 1) {
    	# everything after '|' means alternative-specific coefficients
    	f2 <- formula(paste(f1[[2]], f1[[1]], f.spl[2], sep=''))
    	tms.v <- terms(f2, data = data)
    	tms.v.labels <- attr(tms.v, "term.labels")
    }
    response.var <- deparse(f1[[2]])
    if(is.null(choices)) {
    	choices <- as.vector(unique(data[,response.var]))
    	freq.table <- table(data[,response.var])[choices]
    } else {
    	if(!any(is.element(choices, data[,response.var])))			warning('Values in "choices" do not correspond to values in the "', response.var, '" column of data.')
    	use.idx <- is.element(data[,response.var], choices)
    	freq.table <- table(data[use.idx,response.var])[choices]
    }
    nchoices <- length(choices)
    if(!is.element(base.choice, 1:nchoices)) 
    	stop('base.choice must be an index of choices, i.e. a number between 1 and ', nchoices)
    other.choices.ids <- setdiff(1:nchoices, base.choice)	
    # get short names of varying variables and filter out those in the formula
    varying.vars.long <- colnames(data)[varying]
    varying.vars.used <- c()
    if (length(varying.vars.long) > 0) {
    	pattern <- paste(sep, choices[1], '$', sep='')
    	for (ichoice in 2:nchoices) pattern <- paste(pattern, '|', sep, choices[ichoice], '$', sep='')
    	varying.vars <- gsub(pattern, '', varying.vars.long) 
    	varying.vars.used <- intersect(unique(varying.vars),
    								c(tms.labels, tms.v.labels))
    }
    # get remaining variables
    not.varying <- if(length(varying) > 0) colnames(data)[-varying] else colnames(data)
    not.varying.used <- intersect(not.varying, c(tms.labels, tms.v.labels))
    
    # all variables
    variable.names <- c(not.varying.used, varying.vars.used)
    nvars <- length(variable.names)
    
    c.variables <- tms.labels # not alt.-specific vars
    v.variables <- tms.v.labels # alt.-specific vars
    c.vars.idx <- v.vars.idx <- rep(FALSE, nvars)
    c.vars.idx[is.element(variable.names, c.variables)] <- TRUE    
    v.vars.idx[is.element(variable.names, v.variables)] <- TRUE

    variable_used <- matrix(FALSE, nrow=nchoices, ncol=nvars)    
    variable_used[other.choices.ids,c.vars.idx] <- TRUE
    variable_used[other.choices.ids,v.vars.idx] <- TRUE
    colnames(variable_used) <- variable.names
    rownames(variable_used) <- choices
    
    same.coef <- rep(FALSE, nvars)
    same.coef[c.vars.idx] <- TRUE
    names(same.coef) <- variable.names
    
    connect.var.choices <- function(x) paste(x, choices, sep=sep)
    full.var.names <-  matrix(sapply(variable.names, connect.var.choices), nrow=nchoices, ncol=nvars)
    colnames(full.var.names) <- variable.names
    rownames(full.var.names) <- choices
    
    intercepts <- rep(FALSE, nchoices)
    if(tms.intercept)
    	intercepts[other.choices.ids] <- TRUE
    names(intercepts) <- choices
    	
    return(structure(list(
    					response=response.var, choices=choices, 
    					base.choice=base.choice,
    					variable.used=variable_used,
    					same.coefs=same.coef,
    					full.var.names=full.var.names,
						varying.names=varying.vars.used, 
						intercepts=intercepts, sep=sep,
						frequency=freq.table), class='mnl.spec'))
}


mnl.spec.reduce <- function(spec, var.names=NULL, intercepts=NULL) {
	# Update spec in a way that only variables are kept that are contained in var.names
	# The function assumes that spec is created using all variables and var.names is a subset of those.
	if(!is.null(var.names)) {
		x.names <- var.names
		for (var in colnames(spec$variable.used)) {
			if(!is.element(var, x.names)) {
				if(any(is.element(spec$full.var.names[,var], x.names))) {
					for(choice in spec$choices){
						if(!is.element(spec$full.var.names[choice,var], x.names)) 
							spec$variable.used[choice,var] <- FALSE
					}	
				} else spec$variable.used[,var] <- FALSE
			}	
		}
	}
	if (!is.null(intercepts)) {
		#update intercepts - include only those in 'intercepts'
		spec$intercepts[] <- FALSE
		intercepts <- intercepts[!is.na(intercepts)]
		spec$intercepts[intercepts] <- TRUE
	}
	# remove variables that are not used in any equation
	remove <- (1:dim(spec$variable.used)[2])[colSums(spec$variable.used) == 0]
	if(length(remove)>0) {
		spec$variable.used <- spec$variable.used[,-remove, drop=FALSE]
		spec$full.var.names <- spec$full.var.names[,-remove, drop=FALSE]
		spec$same.coefs <- spec$same.coefs[-remove, drop=FALSE]
	}
	return(spec)
}

summary.mnl.spec <- function(object, ...) {
	cat('\nResponse variable:', object$response)
	cat('\nBase choice name:', object$choices[object$base.choice])
	cat('\nBase choice index:', object$base.choice)
	cat('\nFrequency of alternatives:')
	print(object$frequency)
	cat('\nEquations:')
	cat('\n----------\n')
	vars <- colnames(object$variable.used)
	nvars <- ncol(object$variable.used)
	nchoices <- length(object$choices)
	output <- data.frame(alternative=paste(object$choices, ':', sep=''), intercept=rep('', nchoices),
							stringsAsFactors=FALSE)
	if (nvars > 0) {
		for(ivar in 1:nvars) 
			output <- cbind(output, rep('', nchoices), stringsAsFactors=FALSE)
		colnames(output)[3:(nvars+2)] <- 1:nvars
	}
	for(ichoice in 1:nchoices) {
		choice <- object$choices[ichoice]
		if(object$intercept[ichoice]) output[ichoice,2] <- paste('asc', choice, sep=object$sep)
		if (nvars == 0) next
		if(!any(object$variable.used[ichoice,])) next
		which.vars <- (1:nvars)[object$variable.used[ichoice,]]
		first.var <- min(which.vars)
		for(ivar in which.vars) {
			output[ichoice, ivar+2] <- paste(vars[ivar], 
					if(!object$same.coefs[ivar]) paste(object$sep, choice, sep='') else '', sep='')
			if(ivar == first.var) {
				if(nchar(output[ichoice,2]) > 0) output[ichoice, ivar+2] <- paste('+', output[ichoice, ivar+2])
			} else {
			if(nchar(output[ichoice,ivar+2]) > 0) output[ichoice, ivar+2] <- paste('+', output[ichoice, ivar+2])}
		}

	}
	print(output)
}