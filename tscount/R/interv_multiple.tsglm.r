interv_multiple <- function(...) UseMethod("interv_multiple")

interv_multiple.tsglm <- function(fit, taus=2:length(fit$ts), deltas=c(0,0.8,1), external=FALSE, B=10, signif_level=0.05, start.control_bootstrap, final.control_bootstrap, inter.control_bootstrap, parallel=FALSE, ...){
#Test on (multiple) interventions of UNknown type at a UNknown points in time
  #taus: vector of times which should be considered for the interventions
  #B: Integer (>0). Number of bootstrap samples for estimation of the p-value.
##############################
  #Check and modify arguments:
  tsglm.check(fit)
  if(fit$n_obs != fit$n_eff) stop("Function is not applicable on a maximum likelihood fit which is\nnot based on all observations. Use a fit obtained with argument\n'init.drop=FALSE' instead.")
  ts <- fit$ts
  model <- fit$model
  xreg <- fit$xreg
  n <- length(ts)
  p <- length(model$past_obs)
  q <- length(model$past_mean)
  r <- max(ncol(xreg), 0)
  L <- length(deltas)
  result <- list(interventions=NULL, fit_H0=fit, fit_cleaned=NULL, model_interv=NULL, xreg_interv=NULL, fit_interv=NULL)
  ts_cleaned <- NULL
  track <- list(tau_max=NULL, size=NULL, test_statistic=NULL, p_value=NULL)
  stop_procedure <- FALSE
  iteration <- 0
  repeat{
    iteration <- iteration + 1
    if(iteration > 1) fit <- tsglm(ts=ts, model=model, xreg=xreg, link=fit$link, distr=fit$distr, score=FALSE, info="none", ...) #fit the model to the cleaned time series under the null hypothesis of no (further) intervention from the second iteration on, the model fit for the original time series in the first iteration has been done before
    testresult_single <- vector("list", L)
    for(l in 1:L){
      testresult_single[[l]] <- interv_detect(fit=fit, taus=taus, delta=deltas[l], external=external, B=B, start.control_bootstrap=start.control_bootstrap, final.control_bootstrap=final.control_bootstrap, inter.control_bootstrap=inter.control_bootstrap, parallel=parallel, est_interv=TRUE, ...)   
    }
    test_statistics <- sapply(testresult_single, function(x) x$test_statistic)
    p_values <- sapply(testresult_single, function(x) x$p_value)*L #Bonferroni correction
    taus_max <- sapply(testresult_single, function(x) x$tau_max)
    sizes <- sapply(testresult_single, function(x) x$fit_interv$coefficients[1+p+q+r+1])
    p_values_ordered <- order(-p_values, deltas==1, test_statistics, decreasing=TRUE) #Criteria for choosing the type of intervention: 1. lowest p-value, 2. prefer LS, 3. largest test statistic   
    for(j in 1:L){ #Usually only j=1 (i.e. the intervention with the minimal p-value) is used and then this loop is stopped. Higher j are only used if the intervention with the minimal p-value is significant but the estimated contamination by this intervention is zero at all times (because the estimated size is too small). This is a rare but possible case which would otherwise cause a (possibly endless) re-computation of the p-values with exactly the same time series as before.
      track$tau_max <- rbind(track$tau_max, taus_max)
      track$size <- rbind(track$size, sizes)
      track$test_statistic <- rbind(track$test_statistic, test_statistics)
      track$p_value <- rbind(track$p_value, p_values)
      l_chosen <- p_values_ordered[j]
      if(p_values[l_chosen] > signif_level){
        result$fit_cleaned <- fit
        stop_procedure <- TRUE #this variable will cause an escape from the outer repeat-loop
        break #this escapes from the inner for-loop
      }
      model_interv <- testresult_single[[l_chosen]]$model_interv
      xreg_interv <- testresult_single[[l_chosen]]$xreg_interv
      fit_interv <- testresult_single[[l_chosen]]$fit_interv
      decomposition <- tsglm.interv_decompose(ts=ts, model=model_interv, xreg=xreg_interv, link=fit$link, paramvec=fit_interv$coefficients, isolate=r+1)
      ts <- decomposition$cleaned #use cleaned time series as the time series of the next step              
      result$interventions <- rbind(result$interventions, c(taus_max[l_chosen], deltas[l_chosen], sizes[l_chosen], test_statistics[l_chosen], p_values[l_chosen]))                              
      ts_cleaned <- c(ts_cleaned, list(ts))
      if(all(decomposition$contamination==0)){ #If an intervention is found but its size is too small to have an effect, the cleaned time series is the same as before. In this case a re-computation of the p-values would only lead to different results due to randomness in the bootstrap. Instead the same p-values as before are used and the intervention with the second largest p-value is considered.
        warning(paste("The size of intervention number", nrow(result$interventions), "is too small to have an effect on the time series"))
        if(j==L) stop_procedure <- TRUE #If there are no more possible intervention types left, the whole procedure is stopped.
        next
      }else{ #In the much more likely other case that there is a contamination at at least one time the for-loop is escaped and new p-values are computed.
        break
      }
    }
    if(stop_procedure) break #This escapes from the outer repeat-loop and thus ends the whole procedure if no significant p-value was found or if no more candidates for interventions are left.  
  }
  dimnames(track$tau_max) <- dimnames(track$size) <- dimnames(track$test_statistic) <- dimnames(track$p_value) <- list(1:nrow(track$tau_max), deltas) 
  result$interventions <- rbind(data.frame(tau=integer(0), delta=numeric(0), size=numeric(0), test_statistic=numeric(0), p_value=numeric(0)), result$interventions)
  colnames(result$interventions) <- c("tau", "delta", "size", "test_statistic", "p_value")
  result <- c(result, ts_cleaned=list(ts_cleaned), track=list(track))

  #Fit the model with all interventions which have been found by the procedure (joint estimation of all parameters):
  result$model_interv <- model
  result$xreg_interv <- xreg
  if(nrow(result$interventions) > 0){
    covariate <- interv_covariate(n=n, tau=result$interventions$tau, delta=result$interventions$delta)    
    result$xreg_interv <- cbind(xreg, covariate)
    colnames(result$xreg_interv) <- c(colnames(xreg), colnames(covariate)) 
    result$model_interv$external <- c(model$external, rep(external, nrow(result$interventions)))
    ts <- result$fit_H0$ts #object ts has been replaced by a cleaned time series before and is now set to the original time series
    result$fit_interv <- try(tsglm(ts=ts, model=result$model_interv, xreg=result$xreg_interv, link=fit$link, distr=fit$distr, ...))
  }else{
    result$fit_interv <- result$fit_H0
  }
   
  class(result) <- "interv_multiple"
  return(result)
}
