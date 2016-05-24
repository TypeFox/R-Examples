#' Summary of distance detection function model object
#'
#' Provides a brief summary of data and fitted detection probability model
#' parameters, model selection criterion, and optionally abundance in the
#' covered (sampled) region and its standard error.
#'
#' The argument \code{N} is used to suppress computation of
#' abundance and average detection probability in calls to summarize the
#' \code{ds} and either \code{io.fi} or \code{trial.fi} for summaries of
#' \code{io} and \code{trial} objects respectively which are composed of a
#' \code{ds} model object and a mark-recapture model object. The corresponding
#' print function is called to print the summary results.
#'
#' @export
#' @param object a \code{ddf} model object
#' @param se if TRUE, computes standard errors
#' @param N if TRUE, computes abundance in covered (sampled) region
#' @param fittedmodel full fitted model when called from \code{trial} or
#'   \code{io}
#' @param ddfobj distance sampling object description
#' @param \dots unspecified and unused arguments for S3 consistency
#' @return list of extracted and summarized objects 
#' @note This function is called by the generic function \code{summary} for any
#'   \code{ddf} model object.  Each function can be called directly by the
#'   user, but it is typically safest to use the generic function
#'   \code{summary} which calls the appropriate function based on the type of
#'   \code{ddf} model.
#' @author Jeff Laake
#' @keywords utility
summary.io.fi <- function(object,se=TRUE,N=TRUE,fittedmodel=NULL,
                          ddfobj=NULL,...){
  # Uses: predict.io.fi

  model <- object
  avgp <- function(model,pdot,...){return(pdot)}

  avgp0 <- function(model,pdot,observer,...){
    switch(observer,
           pdot$p1,     # observer = 1
           pdot$p2,     # observer = 2
           pdot$fitted) # observer = 3
  }

  newdat <- model$mr$data
  newdat$distance <- rep(0,length(newdat$distance))

  # get the apex for the gamma function
  if(!is.null(ddfobj) && ddfobj$type=="gamma"){
    newdat$distance <- rep(apex.gamma(ddfobj),2)
  }

  newdat$offsetvalue <- 0
  pred.at0 <- predict(model,newdat)
  if(is.null(fittedmodel)){
    Nhat <- model$Nhat
    pdot <- model$fitted
  }else{
    pdot <- fittedmodel$fitted
    Nhat <- fittedmodel$Nhat
  }

  n <- length(newdat$distance)/2

  timesdetected <- newdat$detected[newdat$observer==1] +
                   newdat$detected[newdat$observer==2]

  n1 <- length(newdat$distance[newdat$observer==1&newdat$detected==1])
  n2 <- length(newdat$distance[newdat$observer==2&newdat$detected==1])
  n3 <- sum(timesdetected==2)

  average.p0.1 <- sum(avgp0(model,pred.at0,observer=1)/pdot)
  average.p0.2 <- sum(avgp0(model,pred.at0,observer=2)/pdot)
  average.p0 <- sum(avgp0(model,pred.at0,observer=3)/pdot)

  ans <- list(n             = n,
              n1            = n1,
              n2            = n2,
              n3            = n3,
              average.p0.1  = average.p0.1/Nhat,
              average.p0.2  = average.p0.2/Nhat,
              average.p0    = average.p0/Nhat,
              cond.det.coef = coef(model),
              aic           = model$criterion)

  # Output estimates relevant to abundance estimation.
  # Note: that these estimates are for the surveyed strip only and are
  # not expanded to some larger region that was sampled.
  if(N){
    ans$average.p <- n/Nhat
    ans$Nhat <- Nhat
  }

  # Compute se for N and p's; the logical N is set when summary.io.fi is being
  # called for method="io.fi". When method="io" the call to summary.io.fi from
  # summary.io sets N=FALSE. So this is used to determine whether the pdots
  # are to be computed from the io.fi model or the io model results. Likewise
  # when method="io", fittedmodel contains the io model object and is not NULL
  if(se){
    if(N | is.null(fittedmodel)){
      se.obj <- calc.se.Np(model, avgp, n, ans$average.p)
      ans$average.p.se <- se.obj$average.p.se
      ans$Nhat.se <- se.obj$Nhat.se
      ans$cv <- se.obj$Nhat.se/Nhat
      cvN <- ans$cv
      Nhatvar.list <- se.obj$Nhatvar.list
      vcov <- se.obj$vcov
    }

    # If there is a nested model, compute the Nhat from the fitted model
    # and use its variance and vector of partials.
    if(!is.null(fittedmodel)){
      Nhat <- fittedmodel$Nhat
      vcov <- solvecov(fittedmodel$hessian)$inv
      Nhatvar.list <- DeltaMethod(fittedmodel$par,NCovered,vcov,.001,
                                  model=fittedmodel,group=TRUE)
      Nhatvar <- Nhatvar.list$variance +
                   sum((1-fittedmodel$fitted)/fittedmodel$fitted^2)
      cvN <- sqrt(Nhatvar)/Nhat
    }

    # Compute se of p_1(0), p_2(0) and p_1(0)+p_2(0)-p_2(0)*p_1(0)
    p.se <- function(model,avgp0,vcov,obs,fittedmodel,Nhatvar.list,this.p,
                     cvN,Nhat,ansthis.p){
      var.pbar.list <- prob.se(model,avgp0,vcov,observer=obs,fittedmodel)
      covar <- t(Nhatvar.list$partial) %*% vcov %*% var.pbar.list$partial +
                var.pbar.list$covar
      var.pbar <- ansthis.p^2 * (cvN^2 + var.pbar.list$var/this.p^2 -
                              2*covar/(this.p*Nhat))
      return(sqrt(var.pbar))
    }

    # this ansthis.p stuff is not okay, need to deal with that
    ans$average.p0.1.se <- p.se(model,avgp0,vcov,1,fittedmodel,Nhatvar.list,
                                average.p0.1,cvN,Nhat,ans$average.p0.1)
    ans$average.p0.2.se <- p.se(model,avgp0,vcov,2,fittedmodel,Nhatvar.list,
                                average.p0.2,cvN,Nhat,ans$average.p0.2)
    ans$average.p0.se   <- p.se(model,avgp0,vcov,3,fittedmodel,Nhatvar.list,
                                average.p0,cvN,Nhat,ans$average.p0)
  }

  ans$mono <- model$ds$aux$mono
  ans$mono.strict <- model$ds$aux$mono.strict

  class(ans) <- "summary.io.fi"
  return(ans)
}
