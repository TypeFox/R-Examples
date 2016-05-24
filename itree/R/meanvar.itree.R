#ALG: itree copies this direclty from rpart.

# SCCS 02/18/97 @(#)meanvar.rpart.s	1.2
meanvar.itree <- function(tree, xlab = "ave(y)", ylab = "ave(deviance)", ...)

{
	if(!inherits(tree, "itree"))
		stop("Not legitimate itree object")
	if(!tree$method=='anova')
		stop("Plot not useful when method is not 'anova'")

	if(!is.null(tree$penalty)){
		warning("meanvar is impacted by using penalties and is NOT comparable to unpenalized meanvar.") 
	}
	
	frame <- tree$frame
	frame <- frame[frame$var == "<leaf>",  ]
	x <- frame$yval
	y <- frame$dev/frame$n
	label <- row.names(frame)
	plot(x, y, xlab = xlab, ylab = ylab, type = "n", ...)
	text(x, y, label)
	invisible(list(x = x, y = y, label = label))
}

meanvar <- function(tree,...) UseMethod('meanvar')
