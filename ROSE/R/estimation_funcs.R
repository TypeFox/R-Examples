#Last modified on 01/28/2014 

##Accuracy measures
accuracy.meas <- function (response, predicted, threshold = 0.5) 
{
	### checks
		if(length(response) != length(predicted)) 
			stop("Response and predicted must have the same length.\n")
		if(length(labels <- levels(factor(response))) != 2) 
			stop("Response must have two levels.\n")
		if(cl <- class(predicted) == "factor" | class(predicted) == "character")
		{
			if(lev <- length(levels(factor(predicted))) != 2)
				stop("predicted must have two levels.\n")
			predicted <- as.numeric(predicted)
		}
	###

	splitted <- split(predicted, response)
	negatives <- splitted[[as.character(labels[1])]]
	n.negatives <- length(negatives)
	positives <- splitted[[as.character(labels[2])]]
	n.positives <- length(positives)

	TP <- sum(positives >= threshold)
	FP <- sum(negatives >= threshold)
	TN <- sum(negatives < threshold)
	FN <- sum(positives < threshold)

	PRECISION <- TP/(TP+FP)
	RECALL <- TP/(TP+FN)
	F <- RECALL*PRECISION/(RECALL+PRECISION)

	out <- list(Call=match.call(), threshold=threshold, precision= PRECISION, recall = RECALL, F=F)
	class(out) <- "accuracy.meas"
	out
}

####print method for accuracy measures
print.accuracy.meas <- function(x, ...)
{
	cat("\n")
	cat("Call: \n")
	print(x$Call)
	cat("\n")
	cat("Examples are labelled as positive when predicted is greater than", x$threshold,"\n")
	cat("\n")
	cat( paste("precision: ", sprintf("%.3f",x$precision),"\n", sep="") )
	cat( paste("recall: ", sprintf("%.3f",x$recall),"\n", sep="") )
	cat( paste("F: ", sprintf("%.3f",x$F),"\n", sep="") )
}

	
######################################################################
##ROC curve and related internal functions
######################################################################

##Roc curve
roc.curve <- function(response, predicted, plotit=TRUE, add.roc=FALSE, n.thresholds=100, ...)
{
	### checks
		if( length(response)!=length(predicted) ) stop("Response and predicted must have the same length.\n")
		if( length( labels <- levels( factor(response) ) ) != 2 ) stop("Response must have two levels.\n")
		if( cl <- class(predicted)=="factor" | class(predicted)=="character" )
		{
			if( lev <- length( levels( factor(predicted) ) ) > 2 ) stop("predicted must have no more than two levels.\n")
			predicted <- as.numeric(factor(predicted))
		}
	###

	thresholds <- sort(unique(predicted))
	ind.thresholds <- round( seq( 1, length(thresholds), len = min(length(thresholds), n.thresholds) ) )
	thresholds <- (c(-Inf, thresholds[ind.thresholds]) + c(thresholds[ind.thresholds], +Inf))*0.5
	splitted <- split(predicted, response)

	negatives <- splitted[[as.character(labels[1])]]
	n.negatives <- length(negatives) 

	positives <- splitted[[as.character(labels[2])]]
	n.positives <- length(positives)

	pts <- sapply(thresholds, f.roc, positives=positives, negatives=negatives, n.positives=n.positives, n.negatives=n.negatives)

	auc <- -sum( ( pts[2,-1] + pts[2,-ncol(pts)] )*diff(pts[1,]) )*0.5
	if(auc<0.5)
	{ 
		auc <- 1-auc 	
		pts[1:2,] <- pts[2:1,]
	}


	if(plotit)
	{
		if(add.roc)
		{
			lines(x=pts[1,], y=pts[2,], ...)
		}
		else
		{
			plot.roc.curve(pts[1,], pts[2,], ...)
			abline(0, 1, col="grey70")
		}
	}

		obj.roc.curve <- list(Call=match.call(), auc=auc, false.positive.rate=pts[1,], true.positive.rate=pts[2,], thresholds=thresholds)
		class(obj.roc.curve) <- "roc.curve"
		obj.roc.curve
}

###print method for roc curve
print.roc.curve <- function(x, ...) 
{
	if( !is.null(x$auc) ) cat( paste("Area under the curve (AUC): ", sprintf("%.3f",x$auc),"\n", sep="") )
}

###summary method for roc curve
summary.roc.curve <- function(object, ...) 
{

	LST <- list( Call=object$Call, auc=object$auc, false.positive.rate=summary(object$false.positive.rate), true.positive.rate=summary(object$true.positive.rate) ) 
	class(LST) <- "summary.roc.curve"
	LST
}

###print method for summary roc curve
print.summary.roc.curve <- function(x, ...) 
{
	cat("\n")
	cat("Call: \n")
	print(x$Call)
	cat("\n")
	cat("Area under the curve (AUC): \n")
	cat(round(x$auc, digits=3),"\n")
	cat("\n")
	cat("False positive rate for evaluated thresholds: \n")
	print(x$false.positive.rate)
	cat("\n")
	cat("True positive rate for evaluated thresholds: \n")
	print(x$true.positive.rate)
	cat("\n")
}

##compute specificity and sensibility
f.roc <- function(x, positives, negatives, n.positives, n.negatives)
{
	c( sum( negatives>x )/n.negatives, sum( positives>=x )/n.positives )
}

###plot the ROC curve with some default parameters in plot()
plot.roc.curve <- function(x, y, ...)
{
	plot.roc.inner(x, y, ...)
}

###plot the ROC curve with some default parameters in plot()
plot.roc.inner <- function(x, y, main="ROC curve", xlab="False positive rate", ylab="True positive rate", xlim=c(0,1), ylim=c(0,1), col=1, type="l", lwd=2, ...)
{
	plot(x, y, main=main, xlab=xlab, ylab=ylab, xlim=xlim, ylim=ylim, col=col, type=type, lwd=lwd,...)
}

