caic.diagnostics <- function(caicObj, which.terms=NULL, which.tests=c("NV","SD","AGE"), plot=TRUE,
                             outlier=3, ultrametric.tol=0.0001, 
                             plot.signif=plot, alpha=0.05, cex.id=0.7,  ...)
{
	
    ## opar <- par(no.readonly=TRUE)
    ## on.exit(par(opar))

	## get the choice of diagnostics
    which.tests <- match.arg(which.tests, c("NV","SD","AGE"), several.ok=TRUE)
   	
	## get the choice of terms
    if(is.null(which.terms)){
        which.terms <- dimnames(attr(terms(formula(caicObj$mod)), "factors"))[[2]] # explanatory terms
     } else {
        if(! all(which.terms %in% dimnames(attr(terms(formula(caicObj$mod)), "factors"))[[2]])) stop("Not all specified terms were present in the model.")
     }
	
	## illogical plot.signif?
	if(! plot ) plot.signif <- FALSE
	
	## get the data
    tab <- caic.table(caicObj, nodalValues=TRUE, validNodes=TRUE, ultrametric.tol)

 	## check whether nodeAge is complete
    if("AGE" %in% which.tests & any(is.na(tab$nodeAge))) {
        warning("Plots of absolute contrasts against node ages requested where node ages are not available.\n",
                "  The tree may not be ultrametric, or the tolerance may need to be adjusted.")
        which.tests <- which.tests[! which.tests=="AGE"]
    }
	
	# if the model already has a finite robust filter then use 
	# that to drop rows before the model fitting
	robust <- attr(caicObj, 'robust')
	nNonrobust <- sum(abs(tab$studentResid) > robust)
	tab <- tab[abs(tab$studentResid) < robust, ]
	
	## test giving pch for outliers
	outlier <- abs(tab$studentResid) >= outlier
    outlier.pch <- ifelse(outlier, 19, 21) 
    outlier.col<- ifelse(outlier, 'red', 'black') 
    
	## borrowed from plot.lm (okay - stolen)
	labels.id <- rownames(tab)
    text.id <- function(x, y, ind, cex=cex.id) {
        labpos <- ifelse(x > mean(range(x)), 2, 4)
        text(x[ind], y[ind], labels.id[ind], xpd = TRUE, 
            pos = labpos, offset = 0.25, cex=cex)
    }

	## get an array to store slopes 
	slopes <- array(NA, dim=c(length(which.tests), 4, length(which.terms)), 
	                dimnames =list(which.tests, c("Estimate", "Std. Error", "t value", "Pr(>|t|)"), which.terms))

	## loop over the terms:
	for(term in which.terms){
		
		## need substitute here (eval(parse(text=paste)))
		ylabExpr <- substitute(expression(abs(plain(' Contrasts in') ~~ VAR ~~ plain(' '))), env=list(VAR=as.name(term)))
		y <- tab[, term]
		 
		## nodal values
		if("NV" %in% which.tests){
			
			x <- tab[, paste("nodal.", term, sep="")]
			
			## model it
	        mod <- lm(y~x)

     		## extract the t-test on the slope
	        slopes['NV',,term] <- summary.lm(mod)$coef[2,]
	
			## plot it and add the line, if requested, if the model is sane, and if it is significant
			if(plot){
			    plot(y ~ x, pch=outlier.pch, col=outlier.col, ylab=ylabExpr, xlab=paste(term, "@ node"), ...)
				if(any(outlier)) text.id(x, y, outlier)
				if(plot.signif && is.finite(slopes['NV',3,term]) && slopes['NV',4,term] <= alpha) abline(mod, col='red')
			}
		}
		
		## expected standard deviation
		if("SD" %in% which.tests){
			
			x <- sqrt(tab$contrVar)
			
			## model it
	        mod <- lm(y~x)

     		## extract the t-test on the slope
	        slopes['SD',,term] <- summary.lm(mod)$coef[2,]
	
			## plot it
			if(plot){
			    plot(y ~ x, pch=outlier.pch, col=outlier.col, ylab=ylabExpr, xlab=paste("SD @ node"), ...)
				if(any(outlier)) text.id(x, y, outlier)
				if(plot.signif && is.finite(slopes['NV',3,term]) && slopes['SD',4,term] <= alpha) abline(mod, col='red')
			}
		}

		## node age
		if("AGE" %in% which.tests){
			
			x <- log(tab$nodeAge)
			
			## model it
	        mod <- lm(y~x)

     		## extract the t-test on the slope
	        slopes['AGE',,term] <- summary.lm(mod)$coef[2,]
	
			## plot it
			if(plot){
			    plot(y ~ x, pch=outlier.pch, col=outlier.col, ylab=ylabExpr, xlab=paste("Age @ node"), ...)
				if(any(outlier)) text.id(x, y, outlier)
				if(plot.signif && is.finite(slopes['NV',3,term]) &&slopes['AGE',4,term] <= alpha) abline(mod, col='red')
			}
		}
	}
    
	# record the number of points not included because of a robustness filter
	attr(slopes, 'robust') <- robust
	attr(slopes, 'nNonrobust') <- nNonrobust
	
    class(slopes) <- "caic.diagnostics"
    return(slopes)
}    
    
print.caic.diagnostics <- function(x, ...){
    
	robust <- attr(x, 'robust')
	nNonrobust<- attr(x, 'nNonrobust')
	cat('\n')
	if(nNonrobust > 0){
		cat("Excluding ", nNonrobust, ifelse(nNonrobust > 1, ' contrasts', ' contrast'), ' with absolute studentised residuals > ', robust, '\n\n', sep='')
	}
	
	
    for(vars in dimnames(x)[[3]]){
		
		## print the coefficients for each - using printCoefmat
		## - this hardcodes the significance stars, so this simply
		##   uses a string of those values
		cat(vars, ':\n')
		
		## - need to present a matrix here
		## - finesse single rows - no easy drop.which
		cf <- x[,,vars]
		if(! is.matrix(cf)) cf <- matrix(cf, ncol=4, dimnames=list(rownames(x), colnames(x)))
        printCoefmat(cf, signif.legend = FALSE, ...)
		cat('\n')
    }

	cat("---\nSignif. codes: 0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1\n")
    invisible(x)
}

caic.robust <- function(caicObj, robust=3){
	
	## Hmm - this is rather tricky because the contrast data is _not_ present
	## in the global environment, which means that the refitting runs foul of
	## looking for 'contrData' (the data frame inside the function)
	
	## ## THIS BEHAVES CREDIBLY BUT THE DATA DETAILS CHANGE IN UNPLEASANT WAYS
	## ## NO IT DOESN'T. CAN'T COPE WITH log(x) - reinterprets 
	## ## I think for the time being this has to be the lm.fit method again.
	## 
	## caicObj$mod$model$robustRes <- abs(rstudent(caicObj$mod)) < outlier
	## robustMod <- update(caicObj$mod, subset=robustRes, data=caicObj$mod$model)
	
	# this doesn't recalculate the studentised residuals - is only relative to the previous model
	
	nonrobust <- abs(caicObj$contrast.data$studentResid) > robust
	valid <- caicObj$contrast.data$validNodes
	contrMD <-  caicObj$contrast.data$contr$explanatory[valid & ! nonrobust,,drop=FALSE]
	contrRS <-  caicObj$contrast.data$contr$response[valid  & ! nonrobust,,drop=FALSE]
	
	attr(contrMD, 'assign') <- attr(caicObj$contrast.data, 'assign') ## replace attributes 
	if(! is.null(attr(caicObj$contrast.data, 'contrasts'))) attr(contrMD, 'contrasts') <- attr(caicObj$contrast.data, 'contrasts')
	
    robustMod <- lm.fit(contrMD, contrRS)
    class(robustMod) <-  "lm"
	robustMod$terms <- attr(caicObj$mod$model, "terms")    
	robustMod$call <- caicObj$mod$call
	
    # put the model.frame in to the lm object so that predict, etc. calls work
    contrData <- as.data.frame(cbind(contrRS, contrMD))
    robustMod$mod$call <- substitute(lm(FORM, data=contrData), list(FORM=formula))
    robustMod$mod$model <- contrData
	
	caicObj$mod <- robustMod
	attr(caicObj, 'robust') <- robust
	
	return(caicObj)
	
}

caic.label <- function(phy, charset=NULL, action="insert", style="CAIC", tips=FALSE){
    
    # OLD2NEW STATUS: CONVERTED...

    if(! inherits(phy, "phylo")) 
         stop("'", deparse(substitute(phy)), "' not of class 'phylo'")
        
    match.arg(action, c("insert", "replace", "append"))
    match.arg(style, c("RLE", "CAIC"))
    
    contrGp <- split(phy$edge[,2], f=phy$edge[,1]) # handily, split retains numeric order not alphabetic...
    caicLab <- character(max(phy$edge)) 
    names(caicLab) <- 1:max(phy$edge)


    if(is.null(charset)) charset <- LETTERS

    # loop the nodes
    for(nd in seq(along=contrGp)){

        parent <- names(contrGp)[nd]
        children <- contrGp[[nd]]
        if(length(children) > length(charset)) stop("Insufficient characters to code polytomies")
        caicLab[children] <- paste(caicLab[parent], charset[1:length(children)], sep="")

    }

    if(style=="RLE"){
        caicLab <- strsplit(caicLab, split="")
        caicLab <- sapply(caicLab, function(X) with(rle(X), paste(ifelse(lengths > 1, lengths, ""), values, sep="", collapse="")))
    }

    # put in the root label
    caicLab[caicLab == ""] <- "@Root"

    # OLD2NEW: intBool <- as.numeric(names(caicLab)) < 0 # internal nodes now from max(phy$edge)-phy$Nnode +1 to  max(phy$edge)
    intBool <- with(phy, 1:max(edge) > (max(edge) - Nnode))

    # insert option changed from an ordered match to edge[,2] in order to preserve the root label
    switch(action, 
        "replace" = { if(tips) phy$tip.label <- caicLab[! intBool]
                      phy$node.label <- caicLab[intBool]},
        "append"  = { if(is.null(phy$node.label)) phy$node.label <- rep("", phy$Nnode)
                      if(tips) phy$tip.label <- paste(phy$tip.label, caicLab[! intBool])
                      phy$node.label <- paste(phy$node.label, caicLab[intBool])},
        "insert"  =   phy$edge.caic.code <- caicLab) #[match(phy$edge[,2], names(caicLab))])

    return(phy)
}

caic.table <- function(caicObj, validNodes=TRUE, nodalValues=FALSE, ultrametric.tol=0.0001, CAIC.codes=FALSE, style="CAIC"){
    
		# simple code to create a data frame of the contrasts from the caic object
		rowID <- names(caicObj$contrast.data$contrVar)
		
        contr <- with(caicObj$contrast.data$contr, cbind(response, explanatory))

        if(nodalValues){
            nv <- with(caicObj$contrast.data$nodalVals, cbind(response, explanatory))
            colnames(nv) <- paste("nodal.", colnames(nv), sep="")
            tab <- as.data.frame(cbind(contr, nv), row.names=rowID)
        } else {
            tab <- as.data.frame(cbind(contr), row.names=rowID)
        }
        
        tab$contrVar <- caicObj$contrast.data$contrVar
        tab$validNodes <- caicObj$contrast.data$validNodes
        tab$nChild <- caicObj$contrast.data$nChild
        tab$nodeDepth <- caicObj$contrast.data$nodeDepth
        if(! is.null(caicObj$data$phy$edge.length) && is.ultrametric(caicObj$data$phy, tol=ultrametric.tol)) {
			tab$nodeAge <- branching.times(caicObj$data$phy) 
		} else { 
			tab$nodeAge <- NA
		} 
		
		# put studentized residuals
        tab$studentResid <- caicObj$contrast.data$studentResid

       if(CAIC.codes){
           Cphy <- caic.label(caicObj$data$phy, style=style)
           tab$CAIC.code <- Cphy$edge.caic.code[match(rowID, names(Cphy$edge.caic.code))]
       }

        if(validNodes) tab <- subset(tab, validNodes, select=-validNodes)
   
       return(tab)

}
