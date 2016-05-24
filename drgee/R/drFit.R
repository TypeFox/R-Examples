drFit <-
    function(object, rootFinder = findRoots, ...){

	if (class(object)!="drgeeData") {
            stop("An object of class \"drgeeData\" is expected")
	}

        if (object$cond & object$olink=="logit") {
            stop("Doubly robust conditional estimation with outcome link
            logit is not possible (yet...).")
        } else if (object$cond) {
            return( dreFitCond(object, omodel = TRUE,
                               rootFinder = rootFinder, ...) )
        } else {
            return( dreFit(object, omodel = TRUE,
                           rootFinder = rootFinder, ...) )
        }
    }
