F.cr.model.matrix <- function( capture, survival, nan, ns ){

    survX <- F.3d.model.matrix( as.formula(survival), nan, ns )
    surv.intercept <- attr(survX, "intercept") == 1
    sur.names <- attr(survX, "variables")
    ny <- length(sur.names)
    
    capX <- F.3d.model.matrix( capture, nan, ns )
    cap.intercept <- attr(capX, "intercept") == 1
    cap.names <- attr(capX, "variables")
    nx <- length(cap.names)

    ans <- list( 
    	capX=capX, 
    	survX=survX, 
    	n.cap.covars=nx, 
    	n.sur.covars=ny, 
    	cap.intercept= cap.intercept, 
    	sur.intercept = surv.intercept, 
    	cap.vars = cap.names, 
    	sur.vars = sur.names)
    ans
}

