#' Average p and average p(0) variance
#'
#' Computes components of variance for average p=n/N and average p(0) with
#' weights based on empirical covariate distribution, if it contains covariates.
#'
#' Need to add equations here as I do not think they exist in any of the texts.
#' These should probably be checked with simulation.
#'
#' @param model ddf model object
#' @param fct function of detection probabilities; currently only
#'   average (over covariates) detection probability p integrated over distance
#'   or average (over covariates) detection probability at distance 0; p(0)
#' @param vcov variance-covariance matrix of parameter estimates
#' @param observer 1,2,3 for primary, secondary, or duplicates for average
#'   p(0); passed to fct
#' @param fittedmodel full fitted ddf model when \code{trial.fi} or
#'   \code{io.fi} is called from \code{trial} or \code{io} respectively
#' @return \item{var}{variance} \item{partial}{partial derivatives of
#'   parameters with respect to fct} \item{covar}{covariance of n and average p
#'   or p(0)}
#' @seealso prob.deriv
#' @author Jeff Laake
prob.se <- function(model,fct,vcov,observer=NULL,fittedmodel=NULL){
  # Functions Used:  DeltaMethod, prob.deriv (in DeltaMethod)

  # First compute variance component due to estimation of detection function
  # parameters. This uses the delta method and produces a v-c matrix if more
  # than one strata
  if(is.null(fittedmodel)){
    vc1.list <- DeltaMethod(model$par,prob.deriv,vcov,.0001,model=model,
                            parfct=fct,observer=observer,fittedmodel=NULL)
  }else{
    vc1.list <- DeltaMethod(fittedmodel$par,prob.deriv,vcov,.0001,model=model,
                           parfct=fct,observer=observer,fittedmodel=fittedmodel)
  }

  vc1 <- vc1.list$variance
  if(!is.null(observer)){
     newdat <- model$mr$data
     newdat$distance <- rep(0,length(newdat$distance))
     newdat$offsetvalue <- 0
     pred.at0 <- predict(model,newdat)
  }else if(!is.null(fittedmodel)){
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

  if(is.null(fittedmodel)){
    pdot <- model$fitted
    if(is.null(observer)){
      vc2 <- sum(fct(model,pdot)^2*(1 - pdot)/pdot^2)
      covar <- sum(fct(model,pdot)*(1 - pdot)/pdot^2)
    }else{
      vc2 <- sum(fct(model,pred.at0,observer)^2*(1 - pdot)/pdot^2)
      covar <- sum(fct(model,pred.at0,observer)*(1 - pdot)/pdot^2)
    }
  }else{
    pdot <- fittedmodel$fitted
    vc2 <- sum(fct(model,pred.at0,observer)^2*(1 - pdot)/pdot^2)
    covar <- sum(fct(model,pred.at0,observer)*(1 - pdot)/pdot^2)
  }

  return(list(var=vc1+vc2,
              partial=vc1.list$partial,
              covar=covar))
}

