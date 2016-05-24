

extractDIC <- function(fit,...){
  UseMethod("extractDIC")
}


extractDIC.merMod <- function(fit,...){
        #REML <- fit@dims["REML"]
#        llik <- logLik(fit, REML)
#        dev <- fit@deviance["ML"]
#        n <- fit@dims["n"]
#        Dhat <- -2 * (llik)
#        pD <- dev - Dhat
#        DIC <- dev + pD[[1]]
#        names(DIC) <- "DIC"
#        return(DIC)      
        is_REML <- isREML(fit)
        llik <- logLik(fit, REML=is_REML)
        dev <- deviance(refitML(fit))
        n <-  getME(fit, "devcomp")$dims["n"]
        Dhat <- -2 * (llik)
        pD <- dev - Dhat
        DIC <- dev + pD[[1]]
        names(DIC) <- "DIC"
        return(DIC)
}




#
#extractAIC.mer <- function(fit,...){
##     REML <- fit@dims["REML"]
##    llik <- logLik(fit, REML)
##    AIC <- AIC(llik)
##    names(AIC) <- "AIC"
##    return(AIC)
#    L <- logLik(refitML(fit))
#    edf <- attr(L,"df")
#    out <- c(edf,-2*L + k*edf)
#    return(out)
#} 
