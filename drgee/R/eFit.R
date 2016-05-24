eFit <-
    function(object, rootFinder = findRoots, ...){

	if (class(object) != "drgeeData") {
            stop("An object of class \"drgeeData\" is expected\n\n")
	}

        if (object$olink == "logit") {

            if (object$elink != "logit") {
                warning("\nAssuming the logit link for the exposure nuisance model\n\n")
            }

            ## Let y and a switch place and replace v with z
            ## and run retrospective logistic regression

            return(oFit(object, inv = TRUE))

            ## If the outcome link is identity or log
        } else {
            if (object$cond) {
                return(dreFitCond(object, omodel = FALSE, rootFinder = rootFinder, ...))
            } else {
                return(dreFit(object, omodel = FALSE, rootFinder = rootFinder, ...))
            }
        }
    }
