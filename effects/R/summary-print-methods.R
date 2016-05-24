# plot, summary, and print methods for effects package
# John Fox and Jangman Hong
#  last modified 2012-11-30 by J. Fox
#  29 June 2011 added grid, rotx and roty arguments to the two plot methods
#   by S. Weisberg
#  21 Dec 2012 modest modification of empty cells with crossed factors
#  2013-01-17: Added factor.ci.style arg to plot.eff() and plot.effpoly(). J. Fox
#  2013-01-18: Added CI bars to multiline plots with factor.ci.style="bars"
#  2013-01-19: Renamed 'factor.ci.style' to 'ci.style'.  Added a 'none' option
#   extended to variate terms if multiline=TRUE, ci.style="bars"
#  2013-01-30: scale arrow "heads" for error bars relative to cex
#  2013-05-31: fixed symbol colors in legends in plot.eff(). J. Fox
#  2013-08-14: fixed bug in restoring warn option. J. Fox
#  2013-08-27: fixed symbols argument for multiline plot in plot.eff(), reported by Ulrike Gromping. J. Fox
#  2013-08-31: fixed handling of ticks.x argument. John
#  2013-09-25: moved plot.eff methods to plot.methods.R for easier work. Michael
#  2013-10-17: added use.splines argument to plot.effpoly.  Sandy


summary.eff <- function(object, type=c("response", "link"), ...){
	result <- list()
	result$header <- paste("\n", gsub(":", "*", object$term), 'effect\n')
	result$offset <- object$offset
	type <- match.arg(type)
	if (type == "response") {
		object$fit <- object$transformation$inverse(object$fit)
		if (!is.null(object$confidence.level)){
			object$lower <- object$transformation$inverse(object$lower)
			object$upper <- object$transformation$inverse(object$upper)
		}
	}
	result$effect <- array(object$fit,     
			dim=sapply(object$variables, function(x) length(x$levels)),
			dimnames=lapply(object$variables, function(x) x$levels))
	if (!is.null(object$se)){
		result$lower.header <- paste('\n Lower', round(100*object$confidence.level, 2), 
				'Percent Confidence Limits\n')
		result$lower <- array(object$lower,   
				dim=sapply(object$variables, function(x) length(x$levels)),
				dimnames=lapply(object$variables, function(x) x$levels))
		result$upper.header <- paste('\n Upper', round(100*object$confidence.level, 2),
				'Percent Confidence Limits\n')
		result$upper <- array(object$upper,   
				dim=sapply(object$variables, function(x) length(x$levels)),
				dimnames=lapply(object$variables, function(x) x$levels))
	}
	if (object$discrepancy > 1e-3) result$warning <- paste("\nWarning: There is an average discrepancy of", 
				round(object$discrepancy, 3),
				"percent \n     in the 'safe' predictions for effect", object$term, '\n')
	class(result) <- "summary.eff"
	result
}

print.summary.eff <- function(x, ...){
	cat(x$header)
	if (x$offset != 0) cat("\noffset = ", x$offset, "\n\n")
	print(x$effect, ...)
	if (!is.null(x$lower)){
		cat(x$lower.header)
		print(x$lower, ...)
		cat(x$upper.header)
		print(x$upper, ...)
	}
	if (!is.null(x$thresholds)){
		cat("\nThresholds:\n")
		print(x$thresholds, ...)
	}
	if (!is.null(x$warning)) cat(x$warning)
	invisible(x)
}

print.eff <- function(x, type=c("response", "link"), ...){
	cat(paste("\n", gsub(":", "*", x$term), 'effect\n'))
	if (x$offset != 0) cat("\noffset = ", x$offset, "\n\n")
	type <- match.arg(type)
	if (type == "response") x$fit <- x$transformation$inverse(x$fit)
	table <- array(x$fit,     
			dim=sapply(x$variables, function(x) length(x$levels)),
			dimnames=lapply(x$variables, function(x) x$levels))
	print(table, ...)
	if (x$discrepancy > 1e-3) cat(paste("\nWarning: There is an average discrepancy of", 
						round(x$discrepancy, 3),
						"percent \n     in the 'safe' predictions for effect", x$term, '\n'))
	invisible(x)
}

print.efflist <- function(x, ...){
	cat(" model: ")
	form <- x[[1]]$formula
	attributes(form) <- NULL
	print(form)
	for (effect in names(x)){
		print(x[[effect]], ...)
	}
	invisible(x) 
}

summary.efflist <- function(object, ...){
	cat(" model: ")
	form <- object[[1]]$formula
	attributes(form) <- NULL
	print(form)
	for (effect in names(object)){
		print(summary(object[[effect]], ...))
	}
	invisible(NULL) 
}


print.effpoly <- function(x, type=c("probability", "logits"), ...){
	type <- match.arg(type)
	x.frame <-as.data.frame(x)
	n.predictors <- length(names(x$x))
	predictors <- names(x.frame)[1:n.predictors]
	y.lev <- x$y.lev
	ylevel.names <- make.names(paste("prob",y.lev))
	colnames(x$prob) <- colnames(x$logit) <- ylevel.names
	y.categories <- matrix(0, nrow=length(x.frame[,predictors[1]]), ncol=length(y.lev))
	for (i in 1:length(y.lev)){
		level <- which(colnames(x$prob)[i] == ylevel.names)
		y.categories[,i] <-  rep(y.lev[level], length(y.categories[,i]))
	}
	y.categories <- as.vector(y.categories)
	y.categories <- factor(y.categories)
	for (i in 1:length(y.lev)){
		cat(paste("\n", gsub(":", "*", x$term), " effect (", type,") for ", y.lev[i], "\n", sep=""))    
		table <- array(if (type == "probability") {x$prob[y.categories==y.lev[i]]}
				else {x$logit[y.categories==y.lev[i]]},     
			dim=sapply(x$variables, function(x) length(x$levels)),
			dimnames=lapply(x$variables, function(x) x$levels))
		print(table, ...)
	}
	if (x$discrepancy > 0.1) cat(paste("\nWarning: There is an average discrepancy of", 
				round(x$discrepancy, 2),
				"percent \n     in the 'safe' predictions for effect", x$term, '\n'))
	invisible(x)
}

summary.effpoly <- function(object, type=c("probability", "logits"), ...){
	type <- match.arg(type)
	x.frame <-as.data.frame(object)
	n.predictors <- length(names(object$x))
	predictors <- names(x.frame)[1:n.predictors]
	y.lev <- object$y.lev
	ylevel.names <- make.names(paste("prob",y.lev))
	colnames(object$prob) <- colnames(object$logit) <- 
		colnames(object$lower.logit) <- colnames(object$upper.logit) <- 
		colnames(object$lower.prob) <- colnames(object$upper.prob)<- ylevel.names
	y.categories <-matrix(0, nrow=length(x.frame[,predictors[1]]), ncol=length(y.lev))
	for (i in 1:length(y.lev)){
		level <- which(colnames(object$prob)[i] == ylevel.names)
		y.categories[,i] <- rep(y.lev[level], length(y.categories[,i]))
	}
	y.categories <- as.vector(y.categories)
	y.categories <- factor(y.categories)
	for (i in 1:length(y.lev)){
		cat(paste("\n", gsub(":", "*", object$term), " effect (" , type, ") for ", y.lev[i], "\n", sep=""))    
		table <- array(if (type == "probability") {object$prob[y.categories==y.lev[i]]}
				else {object$logit[y.categories==y.lev[i]]},     
			dim=sapply(object$variables, function(x) length(x$levels)),
			dimnames=lapply(object$variables, function(x) x$levels))
		print(table, ...)
	}
	if (is.null(object$confidence.level)) return(invisible(NULL))
	for (i in 1:length(y.lev)){
		cat(paste("\n", 'Lower', object$confidence.level*100, 'Percent Confidence Limits for'
				, y.lev[i],'\n'))
		table <- if (type == "probability") object$lower.prob else object$lower.logit
		table <- array(table[y.categories==y.lev[i]],     
			dim=sapply(object$variables, function(x) length(x$levels)),
			dimnames=lapply(object$variables, function(x) x$levels))
		print(table, ...)
	}
	for (i in 1:length(y.lev)){
		cat(paste("\n", 'Upper', object$confidence.level*100, 'Percent Confidence Limits for'
				, y.lev[i],'\n'))
		table <- if (type == "probability") object$upper.prob else object$upper.logit
		table <- array(table[y.categories==y.lev[i]],     
			dim=sapply(object$variables, function(x) length(x$levels)),
			dimnames=lapply(object$variables, function(x) x$levels))
		print(table, ...)
	}
	if (object$discrepancy > 0.1) cat(paste("\nWarning: There is an average discrepancy of", 
				round(object$discrepancy, 2),
				"percent \n     in the 'safe' predictions for effect", object$term, '\n'))
	invisible(NULL)
}

print.efflatent <- function(x, ...){
	cat(paste("\n", gsub(":", "*", x$term), 'effect\n'))
	table <- array(x$fit,     
		dim=sapply(x$variables, function(x) length(x$levels)),
		dimnames=lapply(x$variables, function(x) x$levels))
	print(table, ...)
	cat("\nThresholds:\n")
	print(x$thresholds, ...)
	if (x$discrepancy > 0.1) cat(paste("\nWarning: There is an average discrepancy of", 
				round(x$discrepancy, 3),
				"percent \n     in the 'safe' predictions for effect", x$term, '\n'))
	invisible(x)
}

summary.efflatent <- function(object, ...){
	result <- list()
	result$header <- paste("\n", gsub(":", "*", object$term), 'effect\n')
	result$effect <- array(object$fit,     
		dim=sapply(object$variables, function(x) length(x$levels)),
		dimnames=lapply(object$variables, function(x) x$levels))
	if (!is.null(object$se)){
		result$lower.header <- paste('\n Lower', round(100*object$confidence.level, 2), 
			'Percent Confidence Limits\n')
		result$lower <- array(object$lower,   
			dim=sapply(object$variables, function(x) length(x$levels)),
			dimnames=lapply(object$variables, function(x) x$levels))
		result$upper.header <- paste('\n Upper', round(100*object$confidence.level, 2),
			'Percent Confidence Limits\n')
		result$upper <- array(object$upper,   
			dim=sapply(object$variables, function(x) length(x$levels)),
			dimnames=lapply(object$variables, function(x) x$levels))
	}
	result$thresholds <- object$thresholds
	if (object$discrepancy > 0.1) result$warning <- paste("\nWarning: There is an average discrepancy of", 
			round(object$discrepancy, 3),
			"percent \n     in the 'safe' predictions for effect", object$term, '\n')
	class(result) <- "summary.eff"
	result
}
