#ALG: itree copies this direclty from rpart.

#SCCS  @(#)na.rpart.s	1.5 12/13/99
na.itree <- function(x){
    Terms <- attr(x, 'terms')
    if(!is.null(Terms)) yvar <- attr(Terms, "response") else yvar <- 0L
    if (yvar==0L) {
	xmiss <- is.na(x)
	keep <-  (xmiss %*% rep(1,ncol(xmiss))) < ncol(xmiss)
	}
    else {
	xmiss <- is.na(x[-yvar])
	ymiss <- is.na(x[[yvar]])
	if (is.matrix(ymiss))
	    keep <- ((xmiss %*% rep(1,ncol(xmiss))) < ncol(xmiss)) &
		    ((ymiss %*% rep(1,ncol(ymiss))) == 0 )
	else
	    keep <- ((xmiss %*% rep(1,ncol(xmiss))) < ncol(xmiss)) & !ymiss
	}
    if (all(keep)) x
    else {
	temp <- seq(keep)[!keep]
	names(temp) <- row.names(x)[!keep]
	#the methods for this group are all the same as for na.omit
	class(temp) <- c("na.itree", "omit")
	structure(x[keep,,drop=FALSE], na.action=temp)
	}
}
