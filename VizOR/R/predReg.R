##' Generate an \code{rms:Predict} object or data frame for an ensemble of
##' simulated disease registries
##'
##' Given a function for generating a simulated disease registry, this function
##' generates an ensemble of such registries. It then returns an \code{rms:Predict}
##' object that contains ensemble-averaged predictions and confidence bounds.
##' @title Registry Ensemble Prediction
##' @param genReg A function that returns a simulated registry dataset, taking as
##'               its first parameter the desired size of the simulated registry,
##'               and possibly other parameters passed through via the \code{\dots}
##'               arg
##' @param N Size of generated registries
##' @param M Size of the ensemble
##' @param fit A fitted model usually intended to serve as a template
##'            for a model to be fitted to the simulated registries.
##'            This may be \code{NULL} when \code{do.pred} is provided
##'            explicitly in the call
##' @param adjust.to A list of adjust-to values for the fitted models,
##'                  defaulting to the adjust-to parameters of \code{fit}
##' @param do.pred A function to be run on the data generated during each iteration
##'                of the simulation, generating either an \code{rms:Predict}
##'                object, or a (usually, named) atomic vector
##' @param ... Additional parameters passed to \code{genReg}
##' @return Depending on the return type of \code{do.pred}, either an
##'         \code{rms:Predict} object containing ensemble-averaged predictions
##'         with confidence bounds reflecting their estimated ensemble variance,
##'         or else a data frame collecting the vector returned by \code{do.pred}
##' @author David C. Norris
##' @keywords datagen
##' @export predReg
predReg <- function(genReg, N, M=100, fit=NULL,
                    adjust.to=fit$Design$limits['Adjust to',],
                    do.pred=function(df){
                      fit.call <- fit$call
                      fit.call$data <- quote(df)
                      ## Regrettably, 'Predict' appears to require that
                      ## its argument be in .GlobalEnv
#                      predReg.fit <<- eval(fit.call)
#                      predReg.fit$Design$limits['Adjust to', names(adjust.to)] <<- adjust.to
                      eval(parse(text="predReg.fit <<- eval(fit.call)"))
                      eval(parse(text=
                                 "predReg.fit$Design$limits['Adjust to',names(adjust.to)]<<-adjust.to"))
                      if(eval(parse(text="is(predReg.fit,'lrm')")))
                        eval(parse(text="Predict(predReg.fit, fun=plogis)"))
                      else
                        eval(parse(text="Predict(predReg.fit)"))
                    },
                    ...){
  ## TODO: Consider the design issues highlighted in the discussion below:
  ## In the case of summarizing simulations generating rms:Predict objects at
  ## each iteration, this function performs important, non-trivial calculations
  ## that would be onerous to re-code for each analysis. Can the same be said
  ## for the case of sims generating simple lists or vectors? Maybe not! But
  ## there may well remain some value to implementing a simple, related case
  ## in the same VizOR function. Also, when I refactor this code to encapsulate
  ## most of the parameters in a 'fit' object, this function will at least offer
  ## a convenient unpacking of this data.
  message("Simulating ", M, " cohorts of size ", N)
  ## We obtain a single datadist before generating the ensemble. This avoids
  ## a jagged ensemble which would create problems for summary calculations.
  ## A proportionally modest overhead cost is involved in this, even with size
  ## round(M*N/2).
  ## TODO: In case 'genReg' is a bootstrap-type resampler, the whole population
  ##       should be used here instead. Presently, this intelligence is required
  ##       of the 'genReg' function itself.
  df <- genReg(round(M*N/2), ...)
#  assign(options()$datadist, datadist(df), envir=.GlobalEnv)
  eval(parse(text=paste(options()$datadist, "<<- datadist(df)")))
  ## We also generate in the local scope an appropriately-shaped rms:Predict object,
  ## which is used below as scaffolding on which the returned object is 'rigged up'.
  pred <- do.pred(df)
  predict <- is(pred, "Predict")
  if(!predict && !(is.atomic(pred) && is.vector(pred)))
    stop("Function 'do.pred' should return an rms:Predict or atomic vector result")
  ## To accomodate an intrinsic limitation of 'replicate', we include the '...'
  ## arguments in a named local variable.
  args <- c(list(N), list(...))
  yhat <- yhat2 <- NULL # appease R CMD check
  Ensemble <- replicate(M, { # Generate an ensemble of M simulated studies
    ## By chance, some generated cohorts may choke our model.
    ## Thus, we resample/regenerate until success achieved.
    ## TODO: Count failures and issue a Warning.
    pred <- NULL
    while(is.null(pred)){
      df <- do.call(genReg, args)
      pred <- tryCatch(do.pred(df),
                       error=function(e) return(NULL))
    }
    if(is.null(pred))
      return(NULL)
    if(predict){
      with(pred, cbind(yhat=yhat, yhat2=yhat^2))
    } else {
      pred # Otherwise, simply collect the 'pred' vectors into a matrix
    }
  })
  if(predict){
    ## Rig up a Predict object containing the ensemble-averaged predictions
    ## and confidence bounds reflecting their estimated ensemble variance:
    EnsembleMean <- apply(Ensemble, seq(length(dim(Ensemble))-1), mean) # collapse the last (replication) dim'n
    EnsembleMean <- as.data.frame(EnsembleMean)
    EnsembleMean <- within(EnsembleMean, { # calculate the s.d. and confidence limits yhat
      yhat.sd <- sqrt(yhat2 - yhat^2) # Var(Y) = E(Y^2) - E(Y)^2
      lower <- yhat - 1.96 * yhat.sd
      upper <- yhat + 1.96 * yhat.sd
    })
    pred$yhat <- EnsembleMean$yhat
    pred$upper <- EnsembleMean$upper
    pred$lower <- EnsembleMean$lower
    pred
  } else { # Simply return the collected vectors as a data frame
    as.data.frame(t(Ensemble))
  }
}
