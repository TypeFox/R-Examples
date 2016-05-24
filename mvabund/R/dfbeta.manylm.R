dfbeta.manylm <- dfbeta.manyglm <- 
function (model, infl = manylm.influence(model, do.coef = TRUE), 
    ...) {

    # create a list that contains the matrix of dfbeta for every variable in y
    b <- infl$coefficients
    n.vars <- length(b)
    
    if(is.null(colnames(b[[1]]))) {
	    for(i in 1:n.vars)
		      colnames(b[[i]]) <- variable.names(model)
    }	
    
    if(is.null(rownames(b[[1]]))) {
	     for(i in 1:n.vars)
 		     rownames(b[[i]]) <- rownames(infl$wt.res)
    }

    names(b) <- colnames(model$coefficients)
    
    return(b)

}

