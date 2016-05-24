interv_detect <- function(...) UseMethod("interv_detect")

interv_detect.tsglm <- function(fit, taus=2:length(fit$ts), delta, external=FALSE, B=NULL, info=c("score"), start.control_bootstrap, final.control_bootstrap, inter.control_bootstrap, parallel=FALSE, est_interv=TRUE, ...){
  #Test on a single intervention of known type at a UNknown point in time
  #B: Integer (>0). Number of bootstrap samples for erstimation of the p-value, for B=NULL no p-value is returned.
  #taus: Integer vector (>=0). Times which should be considered for an intervention.
  ##############################
  
  #Check and modify arguments:
  tsglm.check(fit)
  info <- match.arg(info)
  #if(info=="hessian" && fit$link=="link") stop("For a model with logarithmic link argument 'info' needs to be set to \"score\"")
  taus <- sort(unique(taus)) #ensure, that vector of considered time points is sorted and does not include any duplicates
  
  #Function to compute the test statistic:
  compute_test_statistic <- function(model, ts, xreg, link, distr, fit_H0=NULL, taus, delta, external, est_interv=FALSE, stopOnError=FALSE, ...){
    n <- length(ts)
    if(all(ts==0)) return(list(error_message="Time series is constantly zero"))
    if(is.null(fit_H0)){ #only estimate the parameter under the null hypothesis of no intervention if it is not given in the argument fit_H0
      #ML estimation for the model without intervention:
      op <- options(show.error.messages=FALSE) #suppress error messages for the next call
      fit_H0 <- try(tsglm(ts=ts, model=model, xreg=xreg, link=link, distr=distr, score=FALSE, info="none", ...))
      options(op)
      if("try-error" %in% class(fit_H0)) return(list(error_message=fit_H0[[1]]))
    }
    #Add information about intervention effect:
    param_H0_extended <- c(fit_H0$coefficients, 0)
    model_extended <- model
    xreg_extended <- xreg
    xreg_extended <- cbind(xreg, numeric(n))
    model_extended$external <- c(model$external, external)
    condmean_H0 <- tsglm.condmean(link=link, paramvec=param_H0_extended, model=model_extended, ts=ts, xreg=xreg_extended, derivatives=ifelse(info=="hessian", "second", "first"))
    #Compute test statistic for known time for all tau:
    test_statistic_tau <- as.numeric(rep(NA, length(taus)))
    names(test_statistic_tau) <- taus
    for(j in seq(along=taus)){
      xreg_extended <- cbind(xreg, interv_covariate(n=n, tau=taus[j], delta=delta)) #is overwritten for each tau
      loglik <- tsglm.loglik(link=link, paramvec=param_H0_extended, model=model_extended, ts=ts, xreg=xreg_extended, score=TRUE, info=info, condmean=condmean_H0, from=taus[j])
      infomat_corrected <- apply((1/loglik$kappa + fit$sigmasq)*loglik$outerscoreprod, c(2,3), sum)
      test_statistic_temp <- scoretest(Score=loglik$score, G=loglik$info, G1=infomat_corrected, r=1, stopOnError=stopOnError, silent=TRUE)
      if(!is.null(test_statistic_temp$error_message)){
        return(list(error_message=test_statistic_temp$error_message))
      }
      test_statistic_tau[j] <- test_statistic_temp$test_statistic
    }
    index_tau_max <- which.max(test_statistic_tau)
    tau_max <- taus[index_tau_max]
    test_statistic <- test_statistic_tau[index_tau_max]
    covariate <- interv_covariate(n=n, tau=tau_max, delta=delta)
    xreg_extended <- cbind(xreg, covariate) #add intervention with maximum test statistic to the model
    colnames(xreg_extended) <- c(colnames(xreg), colnames(covariate)) 
    result <- list(
      test_statistic=test_statistic,
      test_statistic_tau=test_statistic_tau,
      tau_max=tau_max,
      fit_H0=fit_H0
    )
    if(est_interv){ #ML estimation for the model with intervention at the point in time where the test statistic has its maximum
      fit_interv <- try(tsglm(ts=ts, model=model_extended, xreg=xreg_extended, link=link, distr=distr, score=FALSE, info="none", ...))
      result <- c(result, list(fit_interv=fit_interv))  
    }
    result <- c(result, list(model_interv=model_extended, xreg_interv=xreg_extended)) 
    return(result)
  }
  
  result <- compute_test_statistic(model=fit$model, ts=fit$ts, xreg=fit$xreg, link=fit$link, distr=fit$distr, fit_H0=fit, taus=taus, delta=delta, external=external, est_interv=est_interv, stopOnError=TRUE, ...)
  
  #Bootstrap to compute p-value:
  if(!is.null(B)){
    #Set arguments controlling the estimation in the bootstrap:
    dotdotdot <- list(...)
    if(missing(start.control_bootstrap)){
      if(!is.null(dotdotdot[["start.control"]])) start.control_bootstrap <- dotdotdot[["start.control"]] else start.control_bootstrap <- list()
    }
    bootstrap_noest <- !missing(final.control_bootstrap) && is.null(final.control_bootstrap) #if argument final.control_bootstrap is NULL, then the parameters are not re-estimated for each bootstrap sample but the true parameters used for simulation are used
    if(missing(final.control_bootstrap)){
      if(!is.null(dotdotdot[["start.control"]])) final.control_bootstrap <- dotdotdot[["final.control"]] else final.control_bootstrap <- list()
    }
    if(is.null(final.control_bootstrap)) final.control_bootstrap <- list()      
    if(missing(inter.control_bootstrap)){
      inter.control_bootstrap <- if(!is.null(dotdotdot[["inter.control"]])) dotdotdot[["inter.control"]] else list()
    }
    if(parallel){
      Sapply <- function(X, FUN, ...) parSapply(cl=NULL, X=X, FUN=FUN, ..., simplify=FALSE) 
    }else{
      Sapply <- function(X, FUN, ...) sapply(X=X, FUN=FUN, ..., simplify=FALSE)
    }
    bootstrap <- function(seed=NULL, fit_H0, n, model, xreg, link, distr, taus, delta, external, ...){
      if(!is.null(seed)) set.seed(seed)
      ts.bootstrap <- tsglm.sim(fit=fit_H0)$ts
      fit_H0.bootstrap <- if(bootstrap_noest) fit_H0 else NULL
      dotdotdot <- list(...)
      dotdotdot[names(dotdotdot) %in% c("start.control", "final.control", "inter.control")] <- NULL #remove these arguments to avoid matching multiple arguments in the following call
      result.bootstrap <- do.call(compute_test_statistic, args=c(list(model=model, ts=ts.bootstrap, xreg=xreg, link=link, distr=distr, fit_H0=fit_H0.bootstrap, taus=taus, delta=delta, external=external, est_interv=TRUE, start.control=start.control_bootstrap, final.control=final.control_bootstrap, inter.control=inter.control_bootstrap), dotdotdot))
      result <- ifelse(is.null(result.bootstrap$error_message), result.bootstrap["test_statistic"], result.bootstrap["error_message"])
      return(result)
    }
    bootstrap_test_statistics <- NULL
    bootstrap_errors <- NULL
    B_left <- B
    while(B_left > 0){
      seeds <- sample(1e+8, size=B_left)
      if(B_left==1) Sapply <- function(X, FUN, ...) sapply(X=X, FUN=FUN, ..., simplify=FALSE) #temporary solution for the problem, that parSapply does give an error (Error in cut.default(i, breaks) : 'breaks' are not unique) if applied to a vector of length one (see https://bugs.r-project.org/bugzilla3/show_bug.cgi?id=14898 for a similar error)   
      output.bootstrap <- Sapply(seeds, bootstrap, fit_H0=fit, n=fit$n_obs, model=fit$model, xreg=fit$xreg, link=fit$link, distr=fit$distr, taus=taus, delta=delta, external=external, ...)
      index_errors <- sapply(output.bootstrap, function(x) is.character(x[[1]]))
      bootstrap_test_statistics <- c(bootstrap_test_statistics, unlist(output.bootstrap[!index_errors]))
      B_left <- B - length(bootstrap_test_statistics)
      bootstrap_errors <- c(bootstrap_errors, unlist(output.bootstrap[index_errors]))
      if(length(bootstrap_errors) >= B){
        print(bootstrap_errors)
        stop(paste("Bootstrap has failed more than 'B' times with", B-B_left, "succesful replications and has thus been aborted"))
      } 
    }
    if(length(bootstrap_errors)>0) warning(paste("For", length(bootstrap_errors), "bootstrapped time series no test statistic could be computed\nand new time series was drawn, see error messages in list\nelement '$bootstrap_errors'"))   
    p_value <- sum(bootstrap_test_statistics > result$test_statistic)/(B+1)
    result <- c(result, list(
      p_value=p_value,
      bootstrap_test_statistics=bootstrap_test_statistics,
      bootstrap_errors=bootstrap_errors
    ))
  }
  
  class(result) <- "interv_detect"
  return(result)
}