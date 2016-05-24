# added 2011-11-04 by J. Fox

effects.sem <- function(object, ...) {
	A <- object$A
	m <- object$m
	I <- diag(m)
	endog <- classifyVariables(object$semmod)$endogenous  
	AA <- - A
	diag(AA) <- 1
	Total <- solve(AA) - I
	Indirect <-  Total - A
	result <- list(Total=Total[endog, ], Direct=A[endog, ], Indirect=Indirect[endog, ])
	class(result) <- "semeffects"
	result
}

print.semeffects <- function(x, digits=getOption("digits"), ...){
	cat("\n Total Effects (column on row)\n")
	Total <- x$Total
	Direct <- x$Direct
	Indirect <- x$Indirect
	select <- !(apply(Total, 2, function(x) all( x == 0)) & 
				apply(Direct, 2, function(x) all( x == 0)) & 
				apply(Indirect, 2, function(x) all( x == 0)))
	print(Total[, select], digits=digits)
	cat("\n Direct Effects\n")
	print(Direct[, select], digits=digits)
	cat("\n Indirect Effects\n")
	print(Indirect[, select], digits=digits)
	invisible(x)
}
