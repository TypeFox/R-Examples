
"vcov.geeglm" <- function(object,...){
sobj<-summary(object)
return(sobj$cov.scaled)
}

"modelparm.geeglm" <- function(model, coef.=coef, vcov.=vcov, df = NULL, ...){

    ### extract coefficients and their covariance matrix
    beta <- coef.(object=model)

    sigma <- vcov.(object=model)
    sigma <- as.matrix(sigma)
    if (any(length(beta) != dim(sigma))) { stop("dimensions of coefficients and covariance matrix don't match")}
    ### try to identify non-estimable parameters
    estimable <- rep(TRUE, length(beta))
    if (any(is.na(beta))) {
        estimable[is.na(beta)] <- FALSE
        beta <- beta[estimable]
    }
    if (length(beta) != ncol(sigma) || nrow(sigma) != sum(estimable))
        stop("could not extract coefficients and covariance matrix from ", 
             sQuote("model"))
   if(is.null(df)){df<-model$df.residual}else{if(df<0){stop("df is not positive")}}
    RET <- list(coef = beta, vcov = sigma, df = df, estimable = estimable)
    class(RET) <- "modelparm"
    RET

}
