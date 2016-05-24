## Functions for canonical discriminant HE plots

# version: 0.5-7
# last revised: 10/20/2007 11:07:02 AM
# -- fixed problems in summary.candisc for ndim printed
# -- moved heplot.candisc & heplot3d.candisc to heplot.candisc.R
# -- fixed candisc.mlm and canrsqTable bugs when rank==1
# -- fixed candisc.mlm bug for 1 factor designs
# -- moved plot.candisc to its own file
# -- plot.candisc: added boxplot of canonical scores when ndim==1 
# Now refer to car:::predictor.names for car_2.0
# Now just copy predictor.names from car to avoid namespace problem

## ----------------------------------------------------------------------------
# Provides:
#   candisc(), candisc.mlm() -- calculate candisc object from mlm object
#   print.candisc(), coef.candisc(), plot.candisc() -- candisc methods

## ----------------------------------------------------------------------------
## canonical scores and vectors from an mlm object
## TODO: provide a data.frame / formula method

predictor.names <- function(model, ...) {
	UseMethod("predictor.names")
}

predictor.names.default <- function(model, ...){
	predictors <- attr(terms(model), "variables")
	as.character(predictors[3:length(predictors)])
}

candisc <-
function(mod, ...) UseMethod("candisc")


candisc.mlm <- function(mod, term, type="2", manova, ndim=rank, ...) {
  if (!inherits(mod, "mlm")) stop("Not an mlm object")
	if (missing(manova)) manova <- Anova(mod, type=as.character(type))
	terms <- manova$terms
	if (missing(term)) term <- terms[1]
	E <- manova$SSPE
	H <- manova$SSP[[term]]
	dfe <- manova$error.df
	dfh <- manova$df[[term]]
	Sp <- E / dfe

  # from discproj.R in library(fpc)
  tdecomp <- function(m){
    wm <- eigen(m, symmetric=TRUE)
    p <- ncol(m)
    wmd <- wm$values
    out <- t(wm$vectors %*% diag(sqrt(wmd)))
    out
  }

  Tm <- tdecomp(E)
  eInv <- solve(Tm)
  eHe <- t(eInv) %*% H %*% eInv
  dc <- eigen(eHe, symmetric=TRUE)
  rank <- min(dfh, sum(dc$values>0))
  pct <- 100 * dc$values / sum(dc$values)
  
  if(ndim > rank) {
	  warning(paste("You asked for", ndim, "dimensions, but rank is", rank, ". ndim has been reset to", rank))
	  ndim <- rank
  }

  coeffs.raw <- eInv %*% dc$vectors * sqrt(dfe)
	# should we drop the coeffs corresponding to 0 eigenvalues here or at the end?
  coeffs.raw <- as.matrix(coeffs.raw[,1:ndim])
	rownames(coeffs.raw) <- rownames(H)
	colnames(coeffs.raw) <- cn <- paste('Can', 1:ndim, sep="")

	# These are what SAS calls pooled within-class std. can. coefficients
	coeffs.std <- diag(sqrt(diag(Sp))) %*% coeffs.raw
	rownames(coeffs.std) <- rownames(H)
	colnames(coeffs.std) <- cn

  data <- model.frame(mod)
	Y <- model.response(data)
	Y <- scale(Y, center=TRUE, scale=FALSE)

  scores <- Y %*% coeffs.raw
  scores <- as.matrix(scores[,1:ndim])    
	colnames(scores) <- cn

	# Get the factor(s) corresponding to the term...
  all.factors <- data[, sapply(data, is.factor), drop=FALSE]
	factor.names <- unlist(strsplit(term, ":"))
	factors <- data[factor.names]
  
  # Canonical means for levels of the factor(s)
	means <- aggregate(scores, factors, mean)
  rownames(means) <- do.call(paste,c(means[factor.names],sep=':'))
  means <- means[, -(1:length(factor.names))]
  
	# These are what SAS calls total canonical structure coefficients
	# and what I plot as vectors in my canplot macro
  structure <- cor(Y, scores)

	canrsq <- dc$values[1:ndim] / (1+dc$values[1:ndim])

  # make scores into a data frame containing factors in mod
#	scores <- cbind( model.frame(mod)[,-1], scores )
#### FIXME: scores should also include regressors in the model
#  scores <- cbind( all.factors, as.data.frame(scores) )
  scores <- cbind( model.frame(mod)[predictor.names(mod)], as.data.frame(scores) )
  result <- list(
  	dfh=dfh, dfe=dfe,
  	eigenvalues=dc$values, canrsq=canrsq,
  	pct=pct, rank=rank, ndim=ndim, means=means, 
  	factors=factors, term=term, terms=terms, 
  	coeffs.raw=coeffs.raw, coeffs.std=coeffs.std,
  	structure=structure,
  	scores=scores
  	)
  class(result) <- "candisc"
  result
}

# print method for candisc objects

print.candisc <- function( x, digits=max(getOption("digits") - 2, 3), ...) {
	table <- canrsqTable(x)
    cat(paste("\nCanonical Discriminant Analysis for ", x$term, ":\n\n", sep=""))
    print(table, digits=digits,na.print = "")

	rank <- x$rank
    eigs <- x$eigenvalues[1:rank]
    tests <- seqWilks(eigs, rank, x$dfh, x$dfe)
    tests <- structure(as.data.frame(tests), 
    heading = paste("\nTest of H0: The canonical correlations in the",
                        "\ncurrent row and all that follow are zero\n") , 
        class = c("anova", "data.frame"))
    print(tests)
    invisible(x)      
}


## calculate table of canonical results
canrsqTable <- function( obj ) {
    rank <- obj$rank
    table <- matrix(NA, rank, 5)
    diff <- obj$eigenvalues[1:rank] - c(obj$eigenvalues[2:rank],NA) ## John
    table[,1] <- obj$canrsq
    table[,2] <- obj$eigenvalues[1:rank] ## John
    table[,3] <- ifelse(rank>1, diff, NA)
    table[,4] <- obj$pct[1:rank] ## John
    table[,5] <- cumsum(obj$pct[1:rank]) ## John
    rownames(table) <- 1:rank
    colnames(table) <- c("CanRsq", "Eigenvalue", "Difference", "Percent", "Cumulative")
    table
}

# Args:
#    eig: eigenvalues of HE^{-1}
#    p:  number of response variables
#    df.h:  degrees of freedom for H
#    df.e:  degrees of freedom for E
# See:
#    http://www.gseis.ucla.edu/courses/ed231a1/notes2/can1.html
#    http://www.gseis.ucla.edu/courses/ed231a1/notes3/manova.html

seqWilks <- function (eig, p, df.h, df.e) 
{
    p.full <- length(eig)
    result <- matrix(0, p.full, 4)
    m <- df.e + df.h - ( p.full + df.h + 1)/2
    for (i in seq(p.full)) {
    	test <- prod(1/(1 + eig[i:p.full]))
    	p <- p.full + 1 - i
    	q <- df.h + 1 - i
    	s <- p^2 + q^2 - 5
    	s <- if (s > 0) 
        sqrt(((p * q)^2 - 4)/s)
    	  else 1
    	df1 <- p * q
    	df2 <- m * s - (p * q)/2 + 1
    	result[i,] <- c(test, ((test^(-1/s) - 1) * (df2/df1)), 
        df1, df2)
    }
    result <- cbind(result, pf(result[,2], result[,3], result[,4], lower.tail = FALSE))
    colnames(result) <- c("LR test stat", "approx F", "num Df", "den Df", "Pr(> F)")
    rownames(result) <- 1:p.full
    result
}

## summary method for a candisc object

summary.candisc <- function( object, means=TRUE, scores=FALSE,
	coef=c("std"),   
	ndim,     # default is min(3, rank, #cumsum(object$pct) < 99) unless ndim<rank
  digits=max(getOption("digits") - 2, 4), ...){

	table <- canrsqTable(object)
    cat(paste("\nCanonical Discriminant Analysis for ", object$term, ":\n\n", sep=""))
    print(table, digits=digits,na.print = "")


	if (missing(ndim)) {
		if (object$ndim < object$rank) ndim <- object$ndim    # ndim was set, use that
		else if (object$rank>3)
			ndim <- min ( 3, object$rank, sum(cumsum(object$pct) < 99) )
		else ndim=object$rank
	}
		
    if (means) {
    	cat("\nClass means:\n\n")
		if (ndim<2) print(object$means, digits=digits)
    	else print(as.matrix(object$means[,1:ndim]), digits=digits)
    }
		# allow for printing any of raw, std, structure coeffs, or none (NULL)
		if (!is.null(coef)) {
    	coef = match.arg(coef, c("std", "raw", "structure"), several.ok=TRUE)
    	for (typ in coef) {
    		cat("\n", typ, "coefficients:\n")
    		coeffs = coef( object, type=typ )[,1:ndim]
    		print( coeffs, digits=digits)
    		}
    	}
    
    if (scores) {
    	cat("\nCanonical scores:\n\n")
    	print(object$scores[,paste("Can",1:ndim,sep="")], digits=digits)
    }

    invisible(NULL)
}

## coef method for a candisc object 
coef.candisc <- function( object, type=c("std", "raw", "structure"), ...) {
	  type <- match.arg(type)
	  switch(type,
	  	std = object$coeffs.std,
	  	raw = object$coeffs.raw,
	  	structure = object$structure)
}

