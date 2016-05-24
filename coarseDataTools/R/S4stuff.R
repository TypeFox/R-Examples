##' An S4 Class that stores a fitted coarse data object
##'
##'
##' @description This is the output from \code{dic.fit()}, which contains the important bits of information about the model fit and key options used.
##'
##'
##' @section Slots:
##'  \describe{
##'    \item{\code{ests}:}{Matrix of class \code{"numeric"}. This matrix summarizes the results of fitting the model. Rows correspond to the first parameter , the second parameter and then percentiles specified by the ptiles argument. Columns correspond to the point estimate, the lower and upper bounds on the 95\% confidence interval and the standard error of the point estimate. If the maximization does not converge, this matrix is filled with NAs.}
##'    \item{\code{conv}:}{Object of class \code{"numeric"}. A value of 1 indicates successful convergence; 0 indicates unsuccessful convergence.}
##'    \item{\code{MSG}:}{Object of class \code{"character"}. The error message returned from \code{optim()} if the routine fails to converge.}
##'    \item{\code{loglik}:}{Object of class \code{"numeric"}. Value of the estimated maximum log-likelihood.}
##'    \item{\code{samples}:}{Object of class \code{"data.frame"}. Data frame of bootstrap estimates of parameters (if bootstraps were performed).}
##'    \item{\code{data}:}{Object of class \code{"data.frame"}. Original data used to fit model.}
##'    \item{\code{dist}:}{Object of class \code{"character"}. Failure time distribution fit to data. "L" for log-normal, "G" for gamma, "W" for Weibull, and "E" for Erlang.}
##'    \item{\code{inv.hessian}:}{Object of class \code{"matrix"}. The inverse of the hessian matrix for the likelihood surface at the MLE. Used to determine the standard errors for the percentiles. Note that optimization is done on a transformed scale with all parameters logged for all distributions except the first parameter of the log-normal distribution.}
##'    \item{\code{est.method}:}{Object of class \code{"character"}. Method used for estimation.}
##'    \item{\code{ci.method}:}{Object of class \code{"character"}. Method used for estimation of confidence/credible intervals.}
##'  }
##'
##' @name cd.fit
##' @rdname cd.fit
##' @aliases cd.fit-class
##' @exportClass cd.fit
setClass("cd.fit",
         representation(ests = "matrix",
                        conv = "numeric",
                        MSG="character",
                        loglik="numeric",
                        samples="data.frame",
                        data="data.frame",
                        dist="character",
                        inv.hessian="matrix",
                        est.method="character",
                        ci.method="character"))


##' An S4 Class that stores a MCMC fit coarse data object
##'
##'
##' This is the output from \code{dic.fit.mcmc()}, which contains the important bits of information about the model fit and key options used.
##'
##'
##'@section Slots:
##'  \describe{
##'    \item{\code{ests}:}{Matrix of class \code{"numeric"}. This matrix summarizes the results of fitting the model. Rows correspond to the first parameter , the second paramaeter and then percentiles specified by the ptiles argument. Columns correspond to the point estimate, the lower and upper bounds on the 95\% credible interval and the standard error of the point estimate.}
##'    \item{\code{conv}:}{Object of class \code{"numeric"}. Not used in with \code{dic.fit.mcmc}}
##'    \item{\code{MSG}:}{Object of class \code{"character"}. The error message returned from optim() if the routine fails to converge.}
##'    \item{\code{loglik}:}{Object of class \code{"numeric"}.  Not used in with \code{dic.fit.mcmc}.}
##'    \item{\code{samples}:}{Object of class \code{"data.frame"}. Data frame of posterior draws of parameters.}
##'    \item{\code{data}:}{Object of class \code{"data.frame"}. Original data used to fit model.}
##'    \item{\code{dist}:}{Object of class \code{"character"}. Failure time distribution fit to data. "L" for log-normal, "G" for gamma, "W" for Weibull, and "E" for Erlang.}
##'    \item{\code{inv.hessian}:}{Object of class \code{"matrix"}. Not used in with \code{dic.fit.mcmc}.}
##'    \item{\code{est.method}:}{Object of class \code{"character"}. Method used for estimation.}
##'    \item{\code{ci.method}:}{Object of class \code{"character"}. Method used for estimation of confidence/credible intervals.}
##'  }
##'
##' @name cd.fit.mcmc
##' @rdname cd.fit.mcmc
##' @aliases cd.fit.mcmc-class
##' @exportClass cd.fit.mcmc
setClass("cd.fit.mcmc",
         contains="cd.fit")

## we don't want boots and data printing out all the time
setMethod("show","cd.fit",function(object){
    cat(sprintf("Coarse Data Model Parameter and Quantile Estimates: \n"))
    print(object@ests)
    cat(sprintf("\n-2*Log Likelihood = %.1f \n",-2*object@loglik))
    if (object@conv == 0) warning("This model did not converge. Try different starting values or increase the number of iterations")
    if (object@dist == "L" & nrow(object@inv.hessian) > 0) {
        se.dispersion <- sqrt(object@inv.hessian[2,2] * (exp(object@ests[2,1] + log(object@ests[2,1])))^2) #using delta method
        cat(sprintf("\nNote: dispersion parameter is exp(sdlog). In this case it is %.3f (95%% CI %.3f-%.3f). \n",
                    exp(object@ests[2,1]),
                    exp(object@ests[2,1])-qt(.975, nrow(object@data)-1)*se.dispersion,
                    exp(object@ests[2,1])+qt(.975, nrow(object@data)-1)*se.dispersion
                    ))
    }
})

## we don't want boots and data printing out all the time
setMethod("show","cd.fit.mcmc",function(object){
    cat(sprintf("Coarse Data Model Parameter and Quantile Estimates: \n"))
    print(object@ests)
    if (object@dist == "L") cat(sprintf("\n Note: dispersion parameter is exp(sdlog). In this case it is %.3f. \n",exp(object@ests[2,1])))
    cat(sprintf("Note: please check that the MCMC converged on the target distribution by running multiple chains. MCMC samples are available in the mcmc slot (e.g. my.fit@mcmc) \n"))
})

#' Get the log-likelihood value of a \code{cd.fit} or \code{cd.fit.mcmc} object
#'
#' @param object A \code{cd.fit} or \code{cd.fit.mcmc} object.
#'
#' @return log-likelihood value
#'
#' @rdname logLik-methods
#' @aliases logLik logLik,cd.fit-method
#' @export
#' @importFrom stats logLik
setMethod("logLik",
          "cd.fit",
          function(object){
              object@loglik
          })

##' Plots the estimated survival function with an option to plot the posterior draws or bootstraps behind the fit
##' 
##' @param x \code{cd.fit} or \code{cd.fit.mcmc} object
##' 
##' @param col.main color for plotting the main estimate
##' 
##' @param col.samps color for the samples (should include some alpha transparency)
##' 
##' @param plot.n.samps how many posterior or boostrap samples do you want to plot?
##' 
##' @param add add to exisiting plot?
##' 
##' @param xlim xlim for plot, defaults to the range of the data slot in x
##' 
##' @param ... other options to pass to plot
##' @rdname plot-methods
##' @aliases plot plot,cd.fit-method
##' @importFrom graphics plot
##' @exportMethod plot
setMethod("plot",
          "cd.fit",
          function(x,y,
                   col.main=rgb(230, 85, 13,maxColorValue=255),
                   col.samps=rgb(99,99,99,10,maxColorValue=255),
                   plot.n.samps=200,
                   add=F, xlim=range(x@data),
                   ...) {

              par1 <- x@ests[1,1]
              par2 <- x@ests[2,1]

              ## which samples do we want to plot
              if (nrow(x@samples) > 0){
                  plot.samples <- x@samples[sample(nrow(x@samples),min(plot.n.samps,nrow(x@samples))),]
              }

              xs <-seq(xlim[1],xlim[2],length=100)
              if (!add) plot(-100,-100,xlim=xlim,ylim=c(0,1),...)
              if (x@dist == "L"){
                  if (plot.n.samps> 0 & nrow(x@samples) > 0){
                      apply(plot.samples,1,function(y) lines(xs,plnorm(xs,meanlog=y[1],sdlog=y[2]),col=col.samps))
                  }
                  lines(xs,plnorm(xs,meanlog=par1,sdlog=par2),col=col.main,lwd=2)
              } else if (x@dist == "G" | x@dist == "E"){
                  if (plot.n.samps > 0 & nrow(x@samples) > 0){
                      apply(plot.samples,1,function(y) lines(xs,pgamma(xs,shape=y[1],scale=y[2]),col=col.samps))
                  }
                  lines(xs,pgamma(xs,shape=par1,scale=par2),col=col.main,lwd=2)
              } else if (x@dist == "W"){
                  if (plot.n.samps > 0 & nrow(x@samples) > 0){
                      apply(plot.samples,1,function(y) lines(xs,pweibull(xs,y[1],y[2]),col=col.samps))
                  }
                  lines(xs,pweibull(xs,par1,par2),col=col.main,lwd=2)
              }
          })
