################################################################################
#  Computation of AIC and BIC of many models of class 'parfm'                  #
################################################################################
#                                                                              #
#  The function 'select.parfm' computes the AIC and BIC values                 #
#    of parametric frailty models with                                         #
#    different baseline hazards and                                            #
#    different frailty distributions                                           #
#                                                                              #
#  Its parameters are                                                          #
#   - formula  : a formula object, with the response                           #
#                on the left of a ~ operator,                                  #
#                and the terms on the right.                                   #
#                The response must be a survival object                        #
#                as returned by the Surv function.                             #
#   - cluster  : the name of the variable in data containing cluster IDs       #
#   - strata   : the name of the variable in data containing strata IDs        #
#   - data     : a data.frame in which to interpret the variables named        #
#                in the formula.                                               #
#   - dist     : the vector of the names of the baseline hazards               #
#   - frailty  : the vector of the names of the frailty distribution           #
#   - method   : the optimization method (See optim())                         #
#   - maxit    : the maximum number of iterations (See optim())                #
#   - showtime : show the execution time of each model? (See parfm())          #
#   - correct  : the correction to use in case of many                         #
#                events per cluster to get finite likelihood values.           #
#                When correct!=0 the likelihood is divided by                  #
#                10^(#clusters * correct) for computation,                     #
#                but the value of the log-likelihood in the output             #
#                is the re-adjusted value.                                     #
#                It is used only for models with Positive Stable frailty       #
#                                                                              #
#                                                                              #
#                                                                              #
#  The function returns a list with elements                                   #
#   - AIC : a table with AIC values of the required models                     #
#           with one line   per baseline hazard distribution and               #
#           with one column per frailty         distribution                   #
#   - BIC : a table with BIC values of the required models                     #
#           with one line   per baseline hazard distribution and               #
#           with one column per frailty         distribution                   #
#                                                                              #
#                                                                              #
#   Date: December 21, 2011                                                    #
#   Last modification on: October 17, 2012                                     #
################################################################################

select.parfm <- function(formula,
                         cluster=NULL,
                         strata=NULL,
                         data,
                         inip=NULL,
                         iniFpar=NULL,
                         dist=c("exponential",
                                "weibull",
                                "gompertz",
                                "loglogistic",
                                "lognormal"),
                         frailty=c("none",
                                   "gamma",
                                   "ingau",
                                   "possta",
                                   "lognormal"),
                         method="BFGS",
                         maxit=500,
                         Fparscale=1,
                         correct=0){
  warn <- getOption("warn")
  options(warn=-1)
  
  res <- list(AIC=NULL, BIC=NULL)
  res$AIC <- res$BIC <- matrix(NA, length(dist), length(frailty),
                               dimnames=list(dist, substr(frailty, 1, 6)))
  cat(paste("\n\n### - Parametric frailty models - ###",
            "Progress status:",
            "  'ok' = converged",
            "  'nc' = not converged\n",
            "                Frailty",
            "Baseline           ",
            sep="\n"))
  cat(c(none=" none  ", gamma=" gamma ", ingau=" invGau", 
        possta=" posSta", lognormal=" lognor")[frailty]
      
      )
  for (d in dist) {
    cat("\n")
    cat(c(exponential = "exponential.......",
          weibull     = "Weibull...........",
          gompertz    = "Gompertz..........",
          loglogistic = "loglogistic.......",
          lognormal   = "lognormal.........")[d])
    for (f in frailty) {
      cat("..")
      model <- try(parfm(formula=formula, 
                         cluster=cluster,
                         strata=strata,
                         data=data,
                         inip=inip,
                         iniFpar=iniFpar,
                         dist=d,
                         frailty=f,
                         method=method,
                         maxit=maxit,
                         Fparscale=Fparscale,
                         showtime=FALSE,
                         correct=correct),
                   silent=TRUE)
      if (!("try-error" %in% class(model))){
        res$AIC[d, substr(f, 1, 6)] <- AIC(model)
        res$BIC[d, substr(f, 1, 6)] <- BIC(model)
        cat("ok....")
      } else cat("nc....")
    }
  }
  cat("\n")
  class(res) <- "select.parfm"
  
  options(warn=warn)
  return(res)
}

