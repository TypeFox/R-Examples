#' Distance Sampling Functions
#'
#' Computes values of conditional and unconditional detection functions and
#' probability density functions for for line/point data for single observer or
#' dual observer in any of the 3 configurations (io,trial,rem).
#'
#' Placeholder -- Not functional ----
#'
#' @param model model object
#' @param newdata dataframe at which to compute values; if NULL uses fitting data
#' @param obs 1 or 2 for observer 1 or 2, 3 for duplicates, "." for combined and "All" to return all of the values
#' @param conditional if FALSE, computes p(x) based on distance detection function and if TRUE based on mr detection function
#' @param pdf if FALSE, returns p(x) and if TRUE, returns p(x)*pi(x)/integral p(x)*pi(x)
#' @param finebr fine break values over which line is averaged
#' @return List containing \item{xgrid }{grid of distance values} \item{values }{average detection fct values at the xgrid values}
#' @author Jeff Laake
ds.function <- function(model,newdata=NULL,obs="All",conditional=FALSE,
                        pdf=TRUE,finebr){

  xgrid <- NULL
  linevalues <- NULL

  # Depending on the type of model, setup data to be used for prediction
  # selecting the specified observer and only those detected (detected=1).
  if(model$method=="io"){
    newdat <- model$mr$mr$data
  }else if(model$method=="trial" | model$method=="trial.fi"){
    newdat <- process.data(model$data,model$meta.data)$xmat
    newdat <- newdat[newdat$observer==obs & newdat$detected==1,]
  }else{
    newdat <- model$mr$data
  }

  newdat$offsetvalue=rep(0,dim(newdat)[1])

  # Also, depending on the type of model, get det fct predicted values
  # Note g0 = p(0) and prob.det are integrated values (mu)
  if(model$method=="io" | model$method=="trial"| model$method=="rem"){
    prob.det <- predict(model)$fitted
    newdat$distance <- 0
    ddfobj <- model$ds$ds$aux$ddfobj

    if(ddfobj$type=="gamma"){
      if(model$method=="io"){
        newdat$distance <- rep(apex.gamma(ddfobj),2)
      }else{
        newdat$distance <- as.vector(apex.gamma(ddfobj))
      }
    }

    if(model$method=="trial"){
      g0 <- predict(model$mr$mr,newdat,type="response")
    }else{
      g0 <- predict(model$mr,newdat,integrate=FALSE)$fitted
    }
  }else{
    prob.det <- predict(model,newdat,integrate=TRUE)$fitted
    newdat$distance <- 0
    g0 <- predict(model,newdat,integrate=FALSE)$fitted
  }

  # For each element in the finebr grid, compute a value for x and the
  # appropriate p(x) averaged over all the covariate values
  for (i in 1:(length(finebr)-1)){
    # Compute x as the midpoint of the breaks and store as distance in dataframe
    x <- (finebr[i]+finebr[i+1])/2
    xgrid <- c(xgrid,x)
    newdat$distance <- rep(x,dim(newdat)[1])

    # Based on model compute p(x) from conditional detection function
    if(model$method!="io" & model$method!="rem"){
      cond.det <- predict(model,newdata=newdat,integrate=FALSE)
    }else{
      cond.det <- predict(model$mr,newdata=newdat,integrate=FALSE)
    }

    if(model$method=="io" | model$method=="io.fi"|
       model$method=="rem" | model$method=="rem.fi"){
      p1 <- cond.det$p1
      p2 <- cond.det$p2
    }else{
      p1 <- cond.det$fitted
    }

    # If this is a point independence model (io, trial, rem) compute the
    # delta(x) values; otherwise set to 1
    par <- model$ds$par

    if(model$method=="io" | model$method=="trial" | model$method=="rem"  ){
      detfct.pooled.values <- detfct(newdat$distance[newdat$observer==1],ddfobj,
                               width=model$meta.data$width-model$meta.data$left)
      deltax <- detfct.pooled.values/(cond.det$fitted/g0)
    }else{
      detfct.pooled.values <- cond.det$fitted/g0
      deltax <- rep(1,length(detfct.pooled.values))
    }

    # Depending on observer compute average values for all observed covariate
    # values at specified distance (x)
    if(obs==1){
      linevalues <- c(linevalues,sum(p1*deltax/prob.det)/sum(1/prob.det))
    }else if(obs==2){
      linevalues <- c(linevalues,sum(p2*deltax/prob.det)/sum(1/prob.det))
    }else if(obs==3){
      linevalues <- c(linevalues,
                      sum(g0*detfct.pooled.values/prob.det)/sum(1/prob.det))
    }else{
      linevalues <- c(linevalues,sum(p1*p2*deltax/prob.det)/sum(1/prob.det))
    }
  }

  return(list(xgrid=xgrid,
              values=linevalues))
}
