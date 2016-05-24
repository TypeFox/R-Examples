#' Derivatives for variance of average p and average p(0) variance
#'
#' Used in call to DeltaMethod from prob.se to get first derivatives
#'
#' Need to add equations here as I do not think they exist in any of the texts.
#' These should probably be checked with simulation.
#'
#' @param par detection function parameter values
#' @param model ddf model object
#' @param parfct function of detection probabilities; currently only
#'   average (over covariates) detection probability p integrated over distance
#'   or average (over covariates) detection probability at distance 0; p(0)
#' @param observer 1,2,3 for primary, secondary, or duplicates for average
#'   p(0); passed to fct
#' @param fittedmodel full fitted ddf model when \code{trial.fi} or
#'   \code{io.fi} is called from \code{trial} or \code{io} respectively
#' @return Vector of values from fct at specified parameter values
#' @seealso prob.se
#' @author Jeff Laake
prob.deriv <- function(par,model,parfct,observer=NULL,fittedmodel=NULL){
  # Depending on model method store new parameter values
  # change 1/1/05 jll: used coefficients instead of coef;
  #                    needed for trial.fi and io.fi

  set.par <- function(model,par){
    model$par <- par
    if(model$method!="ds"){
      if(model$method=="io" | model$method=="trial" | model$method=="rem"){
        model$mr$mr$coefficients <-model$par[1:length(model$mr$mr$coefficients)]
        model$ds$par <- model$par[(length(model$mr$mr$coefficients)+1):
                                   length(par)]
        model$ds$ds$aux$ddfobj <-assign.par(model$ds$ds$aux$ddfobj,model$ds$par)
      }else{
        model$mr$coefficients <- model$par
      }
    }
   return(model)
  }

  model <- set.par(model,par)
  if(!is.null(observer)){
    newdat <- model$mr$data
    newdat$distance <- rep(0,length(newdat$distance))
    newdat$offsetvalue <- 0
    pred.at0 <- predict(model,newdat)
  }else{
    if(!is.null(fittedmodel)){
      if(class(fittedmodel)[1]!="rem"){
        newdat <- model$data[model$data$observer==1&model$data$detected==1,]
      }else{
        newdat <- model$data
      }

      newdat <- newdat[newdat$distance <= model$meta.data$width &
                       newdat$distance >= model$meta.data$left,]
      newdat$distance <- rep(0,length(newdat$distance))
      newdat$offsetvalue <- 0
      pred.at0 <- predict(model,newdat)$fitted
    }
  }

  if(is.null(fittedmodel)){
    pdot <- predict(model,compute=TRUE)$fitted
    if(is.null(observer)){
      return(sum(parfct(model,pdot)/pdot))
    }else{
      return(sum(parfct(model,pred.at0,observer)/pdot))
    }
  }else{
    fittedmodel <- set.par(fittedmodel,par)
    pdot <- predict(fittedmodel,compute=TRUE)$fitted

    if(is.null(observer)){
       return(sum(parfct(model,pred.at0)/pdot))
    }else{
       return(sum(parfct(model,pred.at0,observer)/pdot))
    }
  }
}
