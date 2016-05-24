

#******************************************************************************
# with function for objects of class mids
with.mids.1chain <- function(data, expr, ...) {
    # adapted from with.mids.1chain
    call <- match.call()   
    if (class(data) != "mids.1chain"){ 
        stop("The data must have class mids.1chain")     
                                    }
    data <- data$midsobj    
    #-----------------------------
    # original code from with.mids from mice package
    analyses <- as.list(1:data$m)    
    # do the repeated analysis, store the result.
    for (i in 1:data$m) {
        data.i <- mice::complete(data, i)
        analyses[[i]] <- eval( expr = substitute(expr), envir = data.i, 
					enclos = parent.frame())        
        if (is.expression(analyses[[i]])){ 
            analyses[[i]] <- eval(expr = analyses[[i]], 
						envir = data.i, enclos = parent.frame())
										}
    }
    # return the complete data analyses as a list of length nimp
    object <- list(call = call, call1 = data$call, nmis = data$nmis, analyses = analyses)
    oldClass(object) <- c("mira", "matrix")
    #    class(object) <- "mira"
    return(object)
}
#******************************************************************************
