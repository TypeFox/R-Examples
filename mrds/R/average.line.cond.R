#' Average conditional detection function line for plotting
#'
#' For models with covariates the detection probability for each observation
#' can vary.  This function computes an average value for a set of distances to
#' plot an average line to graphically represent the fitted model in plots that
#' compare histograms and the scatter of individual estimated detection
#' probabilities.
#'
#' @param finebr set of fine breaks in distance over which detection function
#'   values are averaged and plotted
#' @param obs value of observer for averaging (1-2 individual observers)
#' @param model ddf model object
#' @return list with 2 elements:
#'   \tabular{ll}{\code{xgrid} \tab vector of gridded distance values \cr
#'   \code{values} \tab vector of average detection function values at
#'    the \code{xgrid} values\cr}
#' @note Internal function called from plot functions for \code{ddf} objects
#' @author Jeff Laake
#' @keywords utility
average.line.cond <- function(finebr,obs,model){
  xgrid <- NULL
  linevalues <- NULL

  # Depending on the type of model, setup data to be used for prediction
  # selecting the specified observer and only those detected (detected=1).
  if(model$method%in%c("io","trial")){
    newdat <- model$mr$mr$data
  }else{
    # if(model$method=="trial" | model$method=="trial.fi"){
    #   newdat=process.data(model$data,model$meta.data)$xmat
    #   newdat=newdat[newdat$observer!=obs & newdat$detected==1,]
    # }else
    if(model$method=="rem.fi"){
      newdat <- model$data
    }else{
      newdat <- model$mr$data
    }
  }

  newdat$offsetvalue <- rep(0,dim(newdat)[1])

  # For each element in the finebr grid, compute a value for x and the
  # appropriate p(x) averaged over all the covariate values
  for (i in 1:(length(finebr)-1)){
    # Compute x as the midpoint of the breaks and store as distance in dataframe
    x <- (finebr[i]+finebr[i+1])/2
    xgrid <- c(xgrid,x)
    newdat$distance <- rep(x,dim(newdat)[1])

    # Based on model compute p(x) from conditional detection function
    if(model$method=="io" | model$method=="trial"|model$method=="rem"){
      cond.det <- predict(model$mr,newdata=newdat,integrate=FALSE)
    }else{
      cond.det <- predict(model,newdata=newdat,integrate=FALSE)
    }

    if(model$method=="io" | model$method=="io.fi"|
       model$method=="rem" | model$method=="rem.fi"){
      p1 <- cond.det$p1
      p2 <- cond.det$p2
    }else{
      p1 <- cond.det$fitted
    }

    # Depending on observer compute average values for all observed
    # covariate values at specified distance (x)
    if(obs==1){
      linevalues <- c(linevalues,mean(p1))
    }else{
      linevalues <- c(linevalues,mean(p2))
    }
  }
  return(list(xgrid=xgrid,
              values=linevalues))
}
