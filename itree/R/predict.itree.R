#ALG: itree mostly copies this direclty from rpart with
# some changes to reflect different column names in the 'frame' object
# for classification problems.

# SCCS @(#)predict.rpart.s	1.11 06/03/01
predict.itree <-
function(object, newdata = list(),
	 type = c("vector", "prob", "class", "matrix"), na.action = na.pass,
         ...)
{ 
    if(!inherits(object, "itree"))  #instead of rpart checking
	stop("Not legitimate itree object")
	
	#ALG: classification or regression?
	is_classification <- (object$method %in% c("class","class_purity","class_extremes"))
	
    mtype <- missing(type)
    type <- match.arg(type)
    if(missing(newdata))
	where <- object$where
    else {
	if(is.null(attr(newdata, "terms"))) {
	    Terms <- delete.response(object$terms)
	    newdata <- model.frame(Terms, newdata, na.action = na.action,
                                      xlev=attr(object, "xlevels"))
            if (!is.null(cl <- attr(Terms, "dataClasses")))
                .checkMFClasses(cl, newdata, TRUE)
        }
	where <- pred.rpart(object, itree.matrix(newdata)) # ALG: changed to itree.matrix
    }
    frame <- object$frame
    ylevels <- attr(object, "ylevels")
    nclass <- length(ylevels)
    if(mtype && nclass > 0L) type <- "prob"
    if(type == "vector" || (type=="matrix" && is.null(frame$yval2))) {
		pred <- frame$yval[where]
		names(pred) <- names(where)
    }
    else if (type == "matrix") {
    	#alg: only get here if yval2 exists... so something with original rpart.
    	pred <- frame$yval2[where,]  
		dimnames(pred) <- list(names(where), NULL)
    }
    else if(type == "class" && nclass > 0L) {
		pred <- factor(ylevels[frame$yval[where]], levels=ylevels)
		names(pred) <- names(where)
    }
    else if (type == "prob" && nclass > 0L) {
    	#alg: changed this. just get node fractions:
    	pred <- frame[where, grep("wt.frac.class",colnames(frame)),drop=FALSE]
		#pred <- frame$yval2[where, 1L + nclass + 1L:nclass, drop = FALSE] #old
		dimnames(pred) <- list(names(where), ylevels)
    }
    else stop("Invalid prediction for itree object")

    # Expand out the missing values in the result
    # But only if operating on the original dataset
    if (missing(newdata) && !is.null(object$na.action))
        pred <- naresid(object$na.action, pred)
    pred
}

