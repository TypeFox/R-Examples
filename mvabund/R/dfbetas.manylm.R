dfbetas.manylm <-
function (model, infl = manylm.influence(model, do.coef = TRUE), 
    ...) {

    xxi 		 <- chol2inv(model$qr$qr, model$qr$rank)
    n.vars 		 <- NCOL(model$coefficients)
    dfbs 		 <- list()
    dfbetai 	 <- dfbeta(model, infl)
    sqrdigxxi	 <- sqrt(diag(xxi))

    for(i in 1:n.vars)
    	dfbs[[i]] <- dfbetai[[i]] /outer(infl$sigma[,i], sqrdigxxi)
    	
    names(dfbs) 	<- colnames(model$coefficients)
    return(dfbs)
}

