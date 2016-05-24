# Plot optimal subsets regressions -- output from regsubsets
# function in leaps package

# last modified 2015-01-27 by J. Fox

subsets <- function(object, ...){
#	if (!require(leaps)) stop("leaps package missing")
	UseMethod("subsets")
}


subsets.regsubsets <- function(object, 
	names=abbreviate(object$xnames, minlength=abbrev), abbrev=1,
	min.size=1, max.size=length(names), legend="interactive",
	statistic=c("bic", "cp", "adjr2", "rsq", "rss"), las=par("las"), cex.subsets=1,
	...) {
	sumry <- summary(object)
	incidence <- sumry$which
	if (object$xnames[1] == "(Intercept)"){
		if (missing(names)) names <- names[-1]
		incidence <- incidence[, -1]
	}
	statistic <- match.arg(statistic)
	stat <- switch(statistic,
		bic = sumry$bic,
		cp = sumry$cp,
		adjr2 = sumry$adjr2,
		rsq = sumry$rsq,
		rss = sumry$rss)
	subset.size <- as.numeric(rownames(incidence))
	select <- subset.size >= min.size & subset.size <= max.size
	subset.size <- subset.size[select]
	stat <- stat[select]
	incidence <- incidence[select, ]
	plot(c(min.size, max.size), range(stat), type="n", xlab="Subset Size", 
		ylab=paste("Statistic:", statistic), las=las, ...)
	for (i in seq(along=stat)){
		adj <- if (subset.size[i] == min.size) 0
			else if (subset.size[i] == max.size) 1
			else .5
		text(subset.size[i], stat[i], 
			do.call("paste", c(as.list(names[incidence[i,]]),sep='-')),
			cex=cex.subsets, adj=adj)
	}
	if (!is.logical(legend)){
		legend(if (!is.na(charmatch(legend[1], "interactive"))) locator(1) 
               else if (is.character(legend)) legend
               else if (is.numeric(legend) && length(legend == 2)) list(x=legend[1], y=legend[2])
               else stop("improper legend argument"),
			legend=apply(cbind(names, names(names)), 1, 
				function(x) do.call("paste", c(as.list(x), sep=": "))), xpd=TRUE)
		return(invisible(NULL))
	}
	else {
		Abbreviation <- names
		return(as.data.frame(Abbreviation))
	}
		
	
}



