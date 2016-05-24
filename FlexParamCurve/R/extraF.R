extraF <-
structure(function # Compare Two \eqn{nlsList} Models Using Extra Sum-of-Squares F-Tests
                   (submodel = 1,
                    ### \eqn{nlsList} model with fewer curve parameters (reduced model)
                    genmodel = 1,
                    ### \eqn{nlsList} model with more curve parameters (general model)
                    warn = TRUE
                    ) {
    ##description<< Function to compare two nested models using extra sum-of-squares F-Tests.                
    ##details<< Models must be entered in the correct order with the reduced model appearing
    ## first in the call and the more general model appearing later. These must be nested models,
    ## i.e. the general model must contain all of the curve parameters in the reduced model and more.
    ## Entering models with the same number of parameters will produce NAs in the output, but
    ## the function will produce seemingly adequate output with non-nested models. The user must
    ## check that models are nested prior to use.
    ##
    ## This function is primarily designed to be called by the model selection functions
    ## \code{\link{pn.modselect.step}} and \code{\link{pn.mod.compare}} but can be used independently.
    ##
    ## Extra sum-of-squares is obtained from: \preformatted{F = (SS1 - SS2)(df1 - df2) / (SS2 / df2)}
    ## where SS = sum-of-squares and df = degrees of freedom, for the more reduced model (1) and the
    ## more general model (2), respectively.
    ##
    ## If the F value is significant then the more general model provides a significant improvement
    ## over the reduced model, but if the models are not significantly different then the reduced
    ## parameter model is to be preferred.
    Envir1 <- try(FPCEnv$env,silent=T)
    env1ck <- try(is.environment(FPCEnv$env),silent=T)
    envck <- try(is.environment(Envir),silent=T)
    env.ck<-2
    if(envck == FALSE | class(envck) == "try-error") env.ck <- (env.ck - 1)
    if(env1ck == FALSE | class(env1ck) == "try-error") env.ck <- (env.ck - 1)
    if(env.ck == 2){
    	if(identical(Envir,Envir1) == F) Envir <- Envir1}
    if(env.ck == 1 & (envck == FALSE | class(envck) == "try-error")) Envir <- Envir1
    FPCEnv$env <- Envir
    FPCEnv$legitmodel <- "legitmodelreset"
    chk <- try(unlist(summary(submodel))["RSE"], silent = TRUE)
    if (class(chk)[1] == "try-error" | class(submodel)[1] != "nlsList" )
        submodel <- 1
    chk1 <- try(unlist(summary(genmodel))["RSE"], silent = TRUE)
    if (class(chk1)[1] == "try-error" | class(genmodel)[1] != "nlsList" )
        genmodel <- 1
    if (class(submodel)[1] == "numeric" | class(genmodel)[1] ==
        "numeric") {
        if (class(submodel)[1] == "numeric" & class(genmodel)[1] ==
            "numeric") {
            FPCEnv$legitmodel <- 1
        output <- data.frame(NA, NA, NA, NA, NA, NA)
        } else {        
         if (class(submodel)[1] == "numeric") {
            is.na(submodel) <- TRUE
            FPCEnv$legitmodel <- genmodel
	 residu <-resid(genmodel)
         origRSE <- as.numeric(unlist(summary(genmodel))["RSE"]) 
	 newRSE <- origRSE  * sqrt( length(residu))/ sqrt(length(residu[!is.na(residu)]))
	 sdf<- unlist(coef(genmodel)[1])
 	 mindfval <- as.numeric(unlist(summary(genmodel)["df"]))
	 mindfval <- round(mean (mindfval[1:length(mindfval)/2], na.rm = TRUE) )
         df1 <- length(residu) - (length(sdf) * mindfval)
         RSSa <- newRSE * df1
         output <- data.frame(NA, df1, NA, NA, RSSa, NA)
         }
         if (class(genmodel)[1] == "numeric") {
            is.na(genmodel) <- TRUE
            FPCEnv$legitmodel <-  submodel
	 residu <-resid(submodel)
         origRSE <- as.numeric(unlist(summary(submodel))["RSE"]) 
	 newRSE <- origRSE  * sqrt( length(residu))/ sqrt(length(residu[!is.na(residu)]))
	 sdf<- unlist(coef(submodel)[1])
 	 mindfval <- as.numeric(unlist(summary(submodel)["df"]))
	 mindfval <- round(mean (mindfval[1:length(mindfval)/2], na.rm = TRUE) )
         df2 <- length(residu) - (length(sdf) * mindfval)
         RSSb <- newRSE * df2
         output <- data.frame(NA, df2, NA, NA, NA, RSSb)    
         }
        }
     } else {
	 residu <-resid(genmodel)
         origRSE <- as.numeric(unlist(summary(genmodel))["RSE"]) 
	 newRSE <- origRSE  * sqrt( length(residu))/ sqrt(length(residu[!is.na(residu)]))
	 sdf<- unlist(coef(genmodel)[1])
 	 mindfval <- as.numeric(unlist(summary(genmodel)["df"]))
	 mindfval <- round(mean (mindfval[1:length(mindfval)/2], na.rm = TRUE) )
         df1 <- length(residu) - (length(sdf) * mindfval)
         RSSa <- newRSE * df1
	 residu <-resid(submodel)
         origRSE <- as.numeric(unlist(summary(submodel))["RSE"]) 
	 newRSE <- origRSE  * sqrt( length(residu))/ sqrt(length(residu[!is.na(residu)]))
	 sdf<- unlist(coef(submodel)[1])
 	 mindfval <- as.numeric(unlist(summary(submodel)["df"]))
	 mindfval <- round(mean (mindfval[1:length(mindfval)/2], na.rm = TRUE) )
         df2 <- length(residu) - (length(sdf) * mindfval)
         RSSb <- newRSE * df2
         F <- (RSSb - RSSa)/(df2 - df1)/(RSSa/df1)
         sigF <- 1 - pf(abs(F), df1, df2)
         output <- data.frame(F, df2 - df1, df2, sigF, RSSa, RSSb)
    }
    names(output) <- c("Fstat", "df_n", "df_d", "P", "RSS_gen",
        "RSS_sub")
    row.names(output) <- paste(paste(paste("test of", as.character(substitute(submodel)),
        sep = " "), "vs.", sep = " "), as.character(substitute(genmodel)),
        sep = " ")
    return(output)
##value<< A \code{\link{data.frame}} listing the names of the models compared, F,
## numerator degrees of freedom,
## demonimator degrees of freedom, P value and the residual sum of squares for both the general
## and reduced models
##references<< Ritz, C. and Streibigg, J. C. (2008) \eqn{Nonlinear regression with R.}
## Springer-Verlag, New York.
##seealso<< \code{\link{extraF.nls}}
##  \code{\link{nlsList}}
##  \code{\link{pn.modselect.step}}
##  \code{\link{pn.mod.compare}}
}
, ex = function(){
   #compare two nested nlsList models (4 vs 8 parameter models)
   modpar(posneg.data$age, posneg.data$mass) #create pnmodelparams for fixed parameters
   # (only first 4 group levels in data used for example's sake)
   subdata<-subset(posneg.data, as.numeric(row.names (posneg.data) ) < 53)
   richardsR1.lis <- nlsList(mass ~ SSposnegRichards(age, Asym = Asym, K = K,
   Infl = Infl, M = M, RAsym = RAsym, Rk = Rk, Ri = Ri, RM = RM, modno = 1)
                        , data = subdata)
   richardsR12.lis <- nlsList(mass ~ SSposnegRichards(age, Asym = Asym, K = K,
   Infl = Infl, M = M, RAsym = 1, Rk = 1, Ri = 1, RM = 1, modno = 12)
                        , data = subdata)
   extraF(richardsR12.lis, richardsR1.lis)
}
)
