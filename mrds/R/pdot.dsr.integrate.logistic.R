#' Compute probability that a object was detected by at least one observer
#'
#' Computes probability that a object was detected by at least one observer
#' (\code{pdot} or p_.) for a logistic detection function that contains
#' distance.
#'
#' @param right either an integration range for binned data (vector of 2) or the rightmost value for integration (from 0 to right)
#' @param width transect width
#' @param beta parameters of logistic detection function
#' @param x data matrix
#' @param integral.numeric set to TRUE unless data are binned (done in this fct) or the model is such that distance is not linear (eg distance^2), If integral.numeric is FALSE it will compute the integral analytically. It should only be FALSE if is.linear.logistic function is TRUE.
#' @param BT FALSE except for the trial configuration; BT stands for Buckland-Turnock who initially proposed a trial configuration for dual observers
#' @param models list of models including \code{g0model}
#' @param GAM Not used at present. The idea was to be able to use a GAM for g(0) portion of detection function; shoudl always be F
#' @param rem only TRUE for the removal configuration but not used and could be removed if pulled from the function calls. Originally thought the pdot integral would differ but it is the same as the io formula. The only thing that differes with removal is that p(2|1)=1. Observer 2 sees everything seen by observer 1,
#' @param point \code{TRUE} for point transects
#'
#' @author Jeff Laake
pdot.dsr.integrate.logistic <- function(right, width, beta, x,
                                        integral.numeric, BT, models, GAM=FALSE,
                                        rem=FALSE, point=FALSE){
  # Functions called:
  # integratelogistic - computes integral (numerically) over x from 0 to width
  #                     of a logistic detection function
  # integratelogisticdup - computes integral (numerically) over x from 0 to
  #                        width of a duplicate logistic detection function
  # integratelogistic.analytic - computes integral (analytically) over x from 0
  #                              to width of a logistic detection function
  # logisticbyx - computes the detection probability with a logistic det fct
  #               for a single set of z covariates but multiple values of x
  # logisticdupbyx - computes the duplicate detection probability with a
  #                  logistic det fct for a single set of z covariates but
  #                  multiple values of x
  # logisticbyz - computes the detection probability at a single x with a
  #               logistic det fct with mulitple sets of z covariates
  #
  # Uniform detection function for g' but g* includes distance
  #
  # If the models are non-linear in distance or data are binned, numerical
  # integation is required for int1 and int2

  # the name right is a bit confusing in that it is either an integration range
  # for binned data or it came be width or left.  In the case that it is width
  # then it is integration from 0 to width. When set to left, the integration
  # is 0 to left and that is subtracted off the integral from 0 to width to
  # yield the integral from left to width. It was done that way to take
  # advantage of the analytical integral function because typically distance is
  # linear in the model.
  if(length(right)>1){
    lower <- right[1]
    right <- right[2]
    integral.numeric <- TRUE
  }else{
    lower <- 0
  }

  if(is.null(x$observer)){
    stop("data must have a column named \"observer\"")
  }

  if(integral.numeric | point){

    if(GAM|length(right)>1){
      is.constant <- FALSE
    }else{
      is.constant <- is.logistic.constant(x[x$observer==1,],
                                          models$g0model,width)
    }
    if(is.constant){
      int1 <- rep(integratelogistic(x=x[x$observer==1,][1,], models, beta,
                                    lower=lower,right, point),
                  nrow(x[x$observer==1,]))
    }else{
      int1 <- rep(NA,nrow(x[x$observer==1,]))
      for(i in 1:nrow(x[x$observer==1,])){
        int1[i] <- integratelogistic(x=(x[x$observer==1,])[i,], models,
                                     beta,lower=lower,right,point)
      }
    }

    if(!BT){
      if(is.logistic.constant(x[x$observer==2,],models$g0model,width)){
        int2 <- rep(integratelogistic(x=x[x$observer==2,][1,], models, beta,
                                      lower=lower,right, point),
                    nrow(x[x$observer==2,]))
      }else{
        int2 <- rep(NA,nrow(x[x$observer==2,]))
        for(i in 1:nrow(x[x$observer==2,])){
          int2[i] <- integratelogistic(x=x[x$observer==2,][i,], models,
                                            beta,lower=lower,right, point)
        }
      }
    }else{
      int2 <- NULL
    }
  }else{
    # If the models are linear in distance, solve int1 and int2 analytically
    int1 <- integratelogistic.analytic(x[x$observer==1,], models=models,
                                       beta=beta, width=right)
    if(!BT){
      int2 <- integratelogistic.analytic(x[x$observer==2,], models=models,
                                         beta=beta, width=right)
    }else{
      int2 <- NULL
    }
  }

  # Numerical integration is always required for int3
  if(!BT){
    if(is.logistic.constant(x[x$observer==1,],models$g0model,width) &
       is.logistic.constant(x[x$observer==2,],models$g0model,width)){

      int3 <- rep(integratelogisticdup(x1=x[x$observer==1,][1,],
                                       x2=x[x$observer==2,][1,],models,beta,
                                       lower=lower,right, point),
                  nrow(x[x$observer==2,]))
    }else{
      int3 <- rep(NA, nrow(x[x$observer==1,]))
      for(i in 1:nrow(x[x$observer==1,])){
        int3[i] <- integratelogisticdup(x1=(x[x$observer==1,])[i,],
                                        x2=(x[x$observer==2,])[i,],
                                        models, beta, lower=lower,
                                        right, point)
      }
    }
    pdot <- int1 + int2 - int3
  }else{
    int3 <- NULL
    pdot <- int1
  }

  if(!point){
    div <- width
  }else{
    div <- width^2
  }

  # Return list of results
  return(list(pdot = pdot/div,
              int1 = int1/div,
              int2 = int2/div,
              int3 = int3/div))
}
