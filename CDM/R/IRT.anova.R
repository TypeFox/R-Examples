
###########################################
# general ANOVA function
IRT.anova <- function( object , ... ){	
    cl <- match.call()
	cl2 <- paste(cl)[-1]
    if (length(list(object, ...)) != 2){ 
        stop("anova method can only be applied for comparison of two models.\n")		
		}
	objects <- list(object, ...)
	model1a <- objects[[1]]
	model2a <- objects[[2]]
	model1 <- IRT.IC( model1a )
	model2 <- IRT.IC( model2a )
	
    dfr1 <- data.frame( "Model" = cl2[1] , 
			"loglike" = model1["loglike"] , 
			"Deviance" = -2*model1["loglike"] )
    dfr1$Npars <- model1["Npars"]
    dfr1$AIC <- model1["AIC"]
    dfr1$BIC <- model1["BIC"]
    dfr2 <- data.frame( "Model" = cl2[2] , 
			"loglike" = model2["loglike"] , 
			"Deviance" = -2*model2["loglike"] )
    dfr2$Npars <- model2["Npars"]
    dfr2$AIC <- model2["AIC"]
    dfr2$BIC <- model2["BIC"]
	
    dfr <- rbind( dfr1 , dfr2 )
    dfr <- dfr[ order( dfr$Npars ), ]
    dfr$Chisq <- NA
    dfr$df <- NA
    dfr$p <- NA
    dfr[1,"Chisq"] <- dfr[1,"Deviance"] - dfr[2,"Deviance"]
    dfr[1,"df"] <- abs( dfr[1,"Npars"] - dfr[2,"Npars"] )
    dfr[ 1, "p" ] <- round( 1 - stats::pchisq( dfr[1,"Chisq"] , df= dfr[1,"df"] ) , 5 )
    for ( vv in 2:( ncol(dfr))){ dfr[,vv] <- round( dfr[,vv] , 5 ) }
	rownames(dfr) <- NULL
    print( dfr )
    invisible(dfr)
	}
#######################################################################	