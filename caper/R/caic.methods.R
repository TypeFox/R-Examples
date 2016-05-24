summary.caic <- function(object, ...){

    summary(object$mod, ...)

}

print.caic <- function(x, ...){

    cat("Phylogenetic Independent Contrasts analysis using:",  attr(x, "contr.method"), ".\n", sep="")
	if(! is.null(attr(x, "macro.method"))) cat("Response values are species rich contrasts using: ", attr(x, "macro.method"), '\n')
    cat("\nPhylogeny: ", attr(x, "phyName"), " (",  length(x$data$phy$tip.label)  ," tips)\n", sep="")
    cat("Data: ",  attr(x, "dataName"), " (",  nrow(x$data$data)  ," rows)\n", sep="")
    cat("Number of valid contrasts: ", sum(x$contrast.data$validNodes), "\n", sep="")
	
	stres <- na.omit(x$contrast.data$studentResid)
	robust <- attr(x, 'robust')
	if(any(abs(stres) > robust)){
		nNonrobust <- sum(abs(stres) > robust)
		cat("Excluding ", nNonrobust, ifelse(nNonrobust > 1, ' contrasts', ' contrast'), ' with absolute studentised residuals > ', robust, '\n', sep='')
	}

    print(summary(x))

}

predict.caic <- function(object, ...){
    
    # need to force the model to get predictions using the contrast table rather than the original data table...
    # don't completely hijack the newdata argument...
    
    dots <- list(...)
    newdataProv <- pmatch(names(dots), "newdata")
    if(all(is.na(newdataProv))) nD <- caic.table(object) else nD <- dots[[newdataProv]]
    predict(object$mod, newdata=nD)
    
    
}

logLik.caic <- function(object, ...){
	
	logLik(object$mod, ...)

}

anova.caic <- function(object, ...){

	## borrowing from anova.lm
	if(length(list(object, ...)) == 1L){
		# no other objects, no other args (scale and test only make sense for multiple models)
		anova(object$mod)
	} else {
		## pass on - having a second function allows the easy interception of test and scale
		## arguments out of the list of objects
		return(anova.caiclist(object, ...))
	}
	
	
}

anova.caiclist <- function(object, ..., scale=0, test='F'){
	
	## ANOVA cares about model types - need to check crunch
	# need to check that the contrast methods are the same
	objects <- list(object, ...)

	objectsClass <- sapply(objects, class)
	if(! all(objectsClass == 'caic')) stop("anova() on mix of 'caic' and non-'caic' objects.")

	objectsContrMethod <- sapply(objects, attr, 'contr.method')
	if(length(unique(objectsContrMethod)) > 1L) stop("anova() on mixed contrast methods")

	objectsMacroMethod <- sapply(objects, attr, 'macro.method')
	if(length(unique(objectsMacroMethod)) > 1L) stop("anova() on mixed macrocaic methods")

	## OK - now pass the mod parts of those object into anova.lmlist()
	mods <- lapply(objects, '[[', 'mod')
	args <- c(mods, list(scale=scale, test=test))
	anv  <- do.call('anova', args)
	# attr(anv, 'heading')[1] <- "Analysis of Variance Table from 'caic' objects\n"
	return(anv)
}

plot.caic <- function(x, ...){
	
	plot(x$mod, ...)
	
}

residuals.caic <- function(object, ...){
	
	residuals(object$mod, ...)
	
}

coef.caic <- function(object, ...){
	
	coef(object$mod, ...)
	
}


# plot.data