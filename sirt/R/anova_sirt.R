
##############################################################
# anova rasch.mml
anova.rasch.mml <- function( object , ... ){
	cl <- base::match.call()
	model1 <- object
	model2 <- list(object , ... )[[2]]
	res <- IRT.anova.sirt( model1 , model2 )
	cl <- base::paste(cl)[-1]
	ind <- base::match( res$Model , c("model1" , "model2") )
	res$Model <- cl[ind]
    print(res)	
	invisible(res)
				}
##############################################################
# anova smirt
anova.smirt <- anova.rasch.mml
##############################################################
# anova rasch.mirtlc
anova.rasch.mirtlc <- anova.rasch.mml
##############################################################
# anova gom
anova.gom <- anova.rasch.mml
##############################################################
# anova rm.facets
anova.rm.facets <- anova.rasch.mml
##############################################################
# anova rm.sdt
anova.rm.sdt <- anova.rasch.mml
##############################################################
# anova prob.guttman
anova.prob.guttman <- anova.rasch.mml
##############################################################



##################################################################
IRT.anova.sirt <- function (object, ...){
    cl <- match.call()
    cl2 <- paste(cl)[-1]
    if (length(list(object, ...)) != 2) {
        stop("anova method can only be applied for comparison of two models.\n")
    }
    objects <- list(object, ...)
    model1a <- objects[[1]]
    model2a <- objects[[2]]
    model1 <- IRT.IC(model1a)
    model2 <- IRT.IC(model2a)
    dfr1 <- base::data.frame(Model = cl2[1], loglike = model1["loglike"], 
        Deviance = -2 * model1["loglike"])
    dfr1$Npars <- model1["Npars"]
    dfr1$AIC <- model1["AIC"]
    dfr1$BIC <- model1["BIC"]
    dfr2 <- data.frame(Model = cl2[2], loglike = model2["loglike"], 
        Deviance = -2 * model2["loglike"])
    dfr2$Npars <- model2["Npars"]
    dfr2$AIC <- model2["AIC"]
    dfr2$BIC <- model2["BIC"]
    dfr <- base::rbind(dfr1, dfr2)
    dfr <- dfr[ base::order(dfr$Npars), ]
    dfr$Chisq <- NA
    dfr$df <- NA
    dfr$p <- NA
    dfr[1, "Chisq"] <- dfr[1, "Deviance"] - dfr[2, "Deviance"]
    dfr[1, "df"] <- abs(dfr[1, "Npars"] - dfr[2, "Npars"])
    dfr[1, "p"] <- round(1 - stats::pchisq(dfr[1, "Chisq"], df = dfr[1, 
        "df"]), 5)
    for (vv in 2:(ncol(dfr))) {
        dfr[, vv] <- round(dfr[, vv], 5)
    }
    rownames(dfr) <- NULL
    # print(dfr)
    invisible(dfr)
}


##############################################################
# Likelihood ratio test for rasch.copula2 objects
anova.rasch.copula2 <- function( object , ... ){
    if (length(list(object, ...)) != 2){ 
        stop("anova method can only be applied for comparison of two models.\n")		
		}
	objects <- list(object, ...)
	model1 <- objects[[1]]
	model2 <- objects[[2]]

	# define some necessary parameters
	model1$AIC <- model1$ic$AIC
	model1$BIC <- model1$ic$BIC
    model1$loglike <- model1$deviance / (-2)
	model1$Npars <- model1$ic$np
	model2$AIC <- model2$ic$AIC
	model2$BIC <- model2$ic$BIC
    model2$loglike <- model2$deviance / (-2)	
	model2$Npars <- model2$ic$np
	# test
    dfr1 <- data.frame( "Model" = "Model 1" , 
		"loglike" = model1$loglike , 
		"Deviance" = -2*model1$loglike )
    dfr1$Npars <- sum(model1$Npars)
    dfr1$AIC <- model1$AIC
    dfr1$BIC <- model1$BIC
    dfr2 <- data.frame( "Model" = "Model 2" , 
		"loglike" = model2$loglike , 	
		"Deviance" = -2*model2$loglike )
    dfr2$Npars <- sum(model2$Npars)
    dfr2$AIC <- model2$AIC
    dfr2$BIC <- model2$BIC
    dfr <- rbind( dfr1 , dfr2 )
    dfr <- dfr[ order( dfr$Npars ), ]
    dfr$Chisq <- NA
    dfr$df <- NA
    dfr$p <- NA
    dfr[1,"Chisq"] <- dfr[1,"Deviance"] - dfr[2,"Deviance"]
    dfr[1,"df"] <- abs( dfr[1,"Npars"] - dfr[2,"Npars"] )
    dfr[ 1, "p" ] <- round( 1 - stats::pchisq( dfr[1,"Chisq"] , df= dfr[1,"df"] ) , 5 )
    for ( vv in 2:( ncol(dfr))){ dfr[,vv] <- round( dfr[,vv] , 5 ) }
    print( dfr )
    invisible(dfr)
            }
anova.rasch.copula3 <- anova.rasch.copula2			
##############################################################