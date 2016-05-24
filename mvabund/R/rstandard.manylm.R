#####################################################################################
# returns the standardized residuals of a multivariate linear model			#
# however standardization is NOT done multivariate, but univariate fashion 		#
#####################################################################################


rstandard.manylm <- function(model, 
                    sd = sqrt(deviance(model)/df.residual(model)), ...) {

    wt.res <- as.matrix(weighted.residuals(model))		    
    hat <- as.vector(diag(model$hat.X))
    n 	<- NROW(wt.res)
    n.vars <- NCOL(wt.res)
    sD 	<- matrix(rep(sd, each=n), nrow=n, ncol=n.vars)
    hatX <- matrix(rep(hat,times=n.vars), nrow=n, ncol=n.vars)
    res <- wt.res /(sD * sqrt(1 - hatX))
    res[is.infinite(res)] <- NaN
    res
}
