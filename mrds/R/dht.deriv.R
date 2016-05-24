#' Computes abundance estimates at specified parameter values using
#' Horvitz-Thompson-like estimator
#'
#' Computes abundance at specified values of parameters for numerical
#' computation of first derivative with respect to parameters in detection
#' function.  An internal function called by DeltaMethod which is invoked by
#' dht.se
#'
#' @param par detection function parameter values
#' @param model ddf model object
#' @param obs observations table
#' @param samples samples table
#' @param options list of options as specified in \code{\link{dht}}
#' @return vector of abundance estimates at values of parameters specified in
#'   par
#' @note Internal function; not intended to be called by user
#' @author Jeff Laake
#' @seealso \code{\link{dht}}, \code{\link{dht.se}}, \code{\link{DeltaMethod}}
#' @keywords utility
dht.deriv <- function(par,model,obs,samples,options=list()){
  # Functions Used:  predict, covered.region.dht, survey.region.dht

  #  Depending on model method store new parameter values
  #  uses coefficients instead of coef; needed for trial.fi and io.fi 
  model$par=par
  if(model$method!="ds"){
    if(model$method%in%c("io","trial","rem")){
      model$mr$mr$coefficients <- model$par[1:length(model$mr$mr$coefficients)]
      model$ds$par <-model$par[(length(model$mr$mr$coefficients)+1):length(par)]
    }else{
      model$mr$coefficients <- model$par
    }
  }

  # Get new predicted detection probabilities and put with observations
  # 1/1/05 jll; change made to fi methods. need integrate=FALSE
  # 1/24/06 jll; must have been on glue when I made the 1/1/05 change; clearly
  # not what I wanted.  You would only not integrate if you didn't want to
  # assume 1/w; otherwise, always want to integrate to get average detection
  # probability
  pdot <- predict(model,integrate=TRUE,compute=TRUE)

  if(!is.null(pdot$fitted)){
    pdot <- pdot$fitted
  }

  # probably should stop here if pdot$fitted is null? anyway...
  # code below would not work if there is a pdot column there already, so remove
  obs$pdot <- NULL

  obs <- merge(obs,data.frame(object=as.numeric(names(model$fitted)),pdot=pdot))

  # Compute covered region abundances by sample depending on value of group
  Nhat.by.sample <- covered.region.dht(obs,samples,options$group)

  # Scale up abundances to survey region
  Nhat.by.sample <- survey.region.dht(Nhat.by.sample, samples,
                                   model$meta.data$width*options$convert.units,
                                   model$meta.data$point)
  Nhat.by.region <- by(Nhat.by.sample$Nhat,Nhat.by.sample$Region.Label,sum)

  # Return vector of predicted abundances

  if(length(Nhat.by.region)>1){
    return(c(as.vector(Nhat.by.region),sum(as.vector(Nhat.by.region))))
  }else{
    return(as.vector(Nhat.by.region))
  }
}
