npmle <- function(data, family = gaussian, maxiter = 500, tol = 1e-4,
                 smooth = TRUE, bass = 0, nmix = NULL)  {
  
  #### For gaussian: estimate=mle, nuisance=std_err
  #### For poisson: estimate=count, nuisance=mean_multiplier 
  #### For binomial: estimate=n.success, nuisance=n.trials
  
  if(class(family)=="newfam") {
    family <- family
  }
  else {
    if (is.character(family))   {
      family <- get(family, mode = "function", envir = parent.frame())
    }
    if (is.function(family))  { 
      family <- family()
    }
    if (is.null(family$family)) {
      stop("family argument is not valid")
    }
  }
  
  switch(family$family,
         gaussian={
            if(sum(data[,2] <= 0) > 0) {
                stop("All reported standard errors must be positive")
            }
            mix.results <- NPestNormal(x=data[,1],std_err=data[,2],maxiter,
                                       tol,nmix) 
         },
         poisson={
            if(sum(data[,1] < 0) > 0) {
                stop("All elements of the first column must be counts")
            }
            if(sum(data[,2] <= 0) > 0) {
                stop("All values of the mean multiplier must be positive")
            }
            mix.results <- NPestPoisson(x = data[,1],eta = data[,2], maxiter, 
                                        tol,nmix)
         },
         binomial={
            if((sum(data[,1] < 0) > 0) | (sum(data[,2] < 0) > 0)) {
                stop("All elements must be nonnegative")
            }
            if(sum(data[,2] < data[,1]) > 0) {
                stop("The number of trials must be greater than or equal
                       to the number of successes")
            }
            mix.results <- NPestBinomial(x = data[,1],ntrials = data[,2], maxiter, 
                                         tol,nmix)
         },
         tdist={
            if(sum(data[,2] <= 0) > 0) {
                stop("All reported standard errors must be positive")
            }
            n <- nrow(data)
            dd <- family$df
            if(length(dd) != 1 & length(dd)!=n) {
                 stop("The number of entered degrees of freedom must be 
                      equal to one or equal to the number of data points.")
            }
            if(length(dd) == 1) {
                dd <- rep(dd,n)
            }
            mix.results <- NPestT(x = data[,1],std_err = data[,2], df = dd,maxiter,
                                  tol,nmix) 
         },
         Gamma={
            stop("gamma family is not available") 
         },
         inverse.gaussian={
            stop("inverse Gaussian is not available")
         },
  )
  o <- order(mix.results$support)
  support <- mix.results$support[o]
  mix.prop <- mix.results$mix.prop[o]
            
  if(smooth)  {
       ### only do the smoothing on the transformed values where prelim.cdf is between 0 and 1
       prelim.cdf <- cumsum(mix.prop)
       oo <- (prelim.cdf > 0) & (prelim.cdf < 1)
       mp <- cumsum(mix.prop[oo])
       ms <- support[oo]
       
       smooth.cdf <- supsmu(ms, qnorm(mp),bass=bass)
       Fhat <- approxfun(smooth.cdf$x,pnorm(smooth.cdf$y),yleft=0,yright=1)
                
       support <- seq(from = min(support),to = max(support),
                      length.out = length(support))
       mix.prop <- c( Fhat(min(support)), diff(Fhat(support)) )
       
       ### due to roundoff error mix.prop sometimes contains small
       ### negative values
       mix.prop[mix.prop < 0] <- 0
       mix.prop <- mix.prop/sum(mix.prop)
       
  }
  else {
       ## No smoothing or interpolation; return the cdf 
       ## as a step function
       Fhat <- stepfun(support,c(0,cumsum(mix.prop)))
  }
  #tmp <- PostProbPois(x = data[,1], eta = data[,2],support,mix.prop)
  #PP <- tmp$postprobs
  ### might do this for each family separately
  
  fhat <- density(support, weights = mix.prop)
  fhat <- approxfun(fhat$x, fhat$y)
  
  ans <- list()
  class(ans) <- "npmix"
  ans$support <- support
  ans$mix.prop <- mix.prop
  ans$loglik <- mix.results$log.lik
  ans$convergence <- mix.results$conv
  ans$numiter <- mix.results$numiter
  ans$Fhat <- Fhat
  ans$fhat <- fhat
  ans$data <- data
  ans$family <- family
  ans$post.mean <- mix.results$post.mean
  #ans$post.prob <- tmp$postprobs
  
  return(ans)
}

plot.npmix <- function(x, density = FALSE, xlim = NULL, ylim = NULL, 
                       main = NULL, xlab = NULL, ylab = NULL, ...)  {
  if(density) {
      ### plot density
      if(is.null(xlim)) {
          xlim <- range(x$support)
      }
      if(is.null(ylim)) {
          ylim <- c(0, 1.03*max(x$fhat(x$support)))
      }
      if(is.null(xlab)) {
          xlab <- "x"
      }
      if(is.null(ylab)) {
          ylab <- "f"
      }
      if(is.null(main)) {
          main <- "smoothed mixture density"
      }
      plot(x$fhat,xlim=xlim,ylim=ylim,xlab=xlab,ylab=ylab,main=main, ...)
  }
  else {
      ### plot CDF
      if(is.null(xlim)) {
          xlim <- range(x$support)
      }
      if(is.null(ylim)) {
          ylim <- c(0,1)
      }
      if(is.null(xlab)) {
          xlab <- "x"
      }
      if(is.null(ylab)) {
          ylab <- "F"
      }
      if(is.null(main)) {
          main <- "Estimated mixture distribution"
      }
      plot(x$Fhat,xlim=xlim,ylim=ylim,xlab=xlab,ylab=ylab,main=main, ...)
   }
}
         
