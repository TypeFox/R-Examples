extraF.nls <-
structure(function # Compare Two \eqn{nls} Models Using Extra Sum-of-Squares F-Tests
                       (submodel,
                        ### \eqn{nls} model with fewer curve parameters (reduced model)
                        genmodel
                        ### \eqn{nls} model with more curve parameters (general model)
                                ) {
    ##description<<Function to compare two nested \code{\link{nls}} models using extra
    ## sum-of-squares F-Tests.                             
    ##details<< Models must be entered in the correct order with the reduced model appearing
    ## first in the call and the more general model appearing later. These must be nested models,
    ## i.e. the general model must contain all of the curve parameters in the reduced model and more.
    ## Entering models with the same number of parameters will produce NAs in the output, but
    ## the function will produce seemingly adequate output with non-nested models. The user must
    ## check that models are nested prior to use.
    ##
    ## This function is not promoted for use in model selection as differences in curves of
    ## different grouping levels in the dataset may be obscured when curves are fitted to the
    ## entire dataset, as in \code{\link{nls}}.
    ##
    ## Extra sum-of-squares is obtained from: \preformatted{F = (SS1 - SS2)(df1 - df2) / (SS2 / df2)}
    ## where SS = sum-of-squares and df = degrees of freedom, for the more reduced model (1) and the
    ## more general model (2), respectively.
    ##
    ## If the F value is significant then the more general model provides a significant improvement
    ## over the reduced model, but if the models are not significantly different then the reduced
    ## parameter model is to be preferred.
    check <- anova(submodel, genmodel)
    df1 <- as.numeric(unlist(summary(genmodel)["df"])[2])
    RSSa <- (as.numeric(unlist(summary(genmodel))["sigma"])^2) *
        df1
    df2 <- as.numeric(unlist(summary(submodel)["df"])[2])
    RSSb <- (as.numeric(unlist(summary(submodel))["sigma"])^2) *
        df2
    F <- (RSSb - RSSa)/(df2 - df1)/(RSSa/df1)
    sigF <- 1 - pf(abs(F), df1, df2)
    output <- data.frame(F, df2 - df1, df2, sigF, RSSa, RSSb)
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
##seealso<< \code{\link{extraF}}
##  \code{\link{nls}}
##  \code{\link{pn.modselect.step}}
##  \code{\link{pn.mod.compare}}
}
, ex = function(){
#fit and compare two nested nls models (7 vs 8 parameter models)
   modpar(posneg.data$age, posneg.data$mass) #create pnmodelparams for fixed parameters
   richardsR1.nls <- nls(mass ~ SSposnegRichards(age, Asym = Asym, K = K,
   Infl = Infl, M = M, RAsym = RAsym, Rk = Rk, Ri = Ri, RM = RM, modno = 1)
                        , data = posneg.data)
   richardsR2.nls <- nls(mass ~ SSposnegRichards(age, Asym = Asym, K = K,
   Infl = Infl, M = M, RAsym = RAsym, Rk = Rk, Ri = Ri, RM = 1, modno = 2)
                        , data = posneg.data)
   extraF.nls(richardsR2.nls, richardsR1.nls)
}
)
