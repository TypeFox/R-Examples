
#############################################################
#                                                           #
#	subset.circular function                            #
#	Author: Claudio Agostinelli                         #
#	E-mail: claudio@unive.it                            #
#	Date: March, 7, 2003                                #
#	Version: 0.1                                        #
#                                                           #
#	Copyright (C) 2003 Claudio Agostinelli              #
#                                                           #
#############################################################
 
### do not work fine yet
subset.circular <- function(x, subset, select, ...) {

	xcircularp <- attr(x, "circularp")
	x <- unclass(x)
	ismatrix <- is.matrix(x)
	x <- as.data.frame(x)
    if(missing(subset))
        r <- TRUE
    else {
        e <- substitute(subset)
        r <- eval(e, x, parent.frame())
        r <- r & !is.na(r)
    }

    if(missing(select))
        vars <- TRUE
    else {
        nl <- as.list(1:ncol(x))
        names(nl) <- names(x)
        vars <- eval(substitute(select),nl, parent.frame())
    }
    x <- x[r,vars,drop=FALSE]

	if (ismatrix) x <- as.matrix(x)

    attr(x, "circularp") <- xcircularp
    attr(x, "class") <- "circular"
    return(x)	
}
