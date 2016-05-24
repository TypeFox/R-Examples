##################################
## FUNCTIONS THAT DO ESTIMATION ##
##################################
##Functions in this file:
## estimate.STmodel         EX:ok
## print.estimateSTmodel    EX:ok
## coef.estimateSTmodel     EX:ok
## estimate                 EX:NO, S3 compliance methods

##' Estimates parameters of the spatio-temporal model using maximum-likelihood,
##' profile maximum likelihood or restricted maximum likelihood (REML). The
##' function uses the \emph{L-BFGS-B} method in \code{\link[stats:optim]{optim}}
##' to maximise \code{\link{loglikeST}}.
##' 
##' The starting point(s) for the optimisation can either contain both
##' regression parameters and log-covariances parameters for a total of
##' \code{loglikeSTdim(object)$nparam} parameters or only contain the
##' log-covariances covariances parameters \cr
##' i.e. \code{loglikeSTdim(object)$nparam.cov} parameters. \cr
##' If regression parameters are given but not needed (\code{type!="f"}) they
##' are dropped; if they are needed but not given they are inferred through a
##' generalised least squares (GLS) computation, obtained by calling
##' \code{\link{predict.STmodel}}.
##' 
##' If multiple starting points are used this function returns all optimisation
##' results, along with an indication of the best result. The best result is
##' determined by first evaluating which of the optimisations have converged.
##' Convergence is determined by checking that the output from
##' \code{\link[stats:optim]{optim}} has \code{convergence==0} and that the
##' \code{hessian} is negative definite, \cr
##' i.e. \code{all(eigen(hessian)$value < -1e-10)}. \cr
##' Among the converged optimisations the one with the highest log-likelihood
##' value is then selected as the best result.
##' 
##' If none of the optimisations have converged the result with the highest
##' log-likelihood value is selected as the best result.
##' 
##' Most of the elements in \code{res.best} (and in \code{res.all[[i]]}) are
##' obtained from \code{\link[stats:optim]{optim}}. The following is a
##' brief description:
##' \describe{
##'   \item{par}{The best set of parameters found.}
##'   \item{value}{Log-likelihood value corresponding to \code{par}.}
##'   \item{counts}{The number of function/gradient calls.}
##'   \item{convergence}{\code{0} indicates successful convergence,
##'                      see \code{\link[stats:optim]{optim}}.}
##'   \item{message}{Additional information returned by
##'                  \code{\link[stats:optim]{optim}}.}
##'   \item{hessian}{A symmetric matrix giving the finite difference Hessian of
##'                  the function \code{par}.}
##'   \item{conv}{A logical variable indicating convergence; \code{TRUE} if
##'               \code{convergence==0} and \code{hessian} is negative definite,
##'               see \code{details} above.}
##'   \item{par.init}{The initial parameters used for this optimisation.}
##'   \item{par.all}{All parameters (both regression and \emph{log}-covariance).
##'                  Identical to \code{par} if \code{type="f"}.}
##'    \item{hessian.all}{The hessian for all parameters (both regression and
##'                       \emph{log}-covariance).\cr
##'                       \strong{NOTE:} Due to computational considerations
##'                       \code{hessian.all} is computed \emph{only} for \cr
##'                       \code{res.best}.}
##' }
##'
##' @title Estimation of the Spatio-Temporal Model
##' 
##' @param object \code{STmodel} object for which to estimate parameters.
##' @param x Vector or matrix of starting point(s) for the optimisation. A
##'   vector will be treated as a single starting point. If \code{x} is a matrix
##'   the optimisation will be run using each column as a separate starting point.
##'   If \code{x} is a single integer then multiple starting points will be
##'   created as a set of constant vectors with the values of each starting point
##'   taken as \code{seq(-5, 5, length.out=x)}. See \code{details} below.
##' @param x.fixed Vector with parameter to be held fixed; parameters marked as
##'   \code{NA} will still be estimated.
##' @param type A single character indicating the type of log-likelihood to use.
##'   Valid options are "f", "p", and "r", for \emph{full}, \emph{profile} or
##'   \emph{restricted maximum likelihood} (REML).
##' @param h,diff.type Step length and type of finite difference to use when
##'   computing gradients, see \code{\link{loglikeSTGrad}}.
##' @param hessian.all If \code{type!="f"} computes hessian (and uncertainties)
##'   for both regression and \emph{log}-covariance parameters, not only for
##'   \emph{log}-covariance parameters. See \code{value} below.
##' @param lower,upper,method Parameter bound and optimisation method,
##'   passed to \code{\link[stats:optim]{optim}}.
##' @param control A list of control parameters for the optimisation.
##'   See \code{\link[stats:optim]{optim}} for details; setting \code{trace}=0
##'   eliminates all ouput.
##' @param restart Number of times to restart each optimisation if
##'   \code{\link[stats:optim]{optim}} fails to converge; can
##'   sometimes resolve issues with L-BFGS-B line search.
##' @param ... Ignored additional arguments.
##' 
##' @return \code{estimateSTmodel} object containing:
##'   \item{res.best}{A list containing the best optimisation result; elements
##'                   are described below. Selection of the best result is described
##'                   in \code{details} above.}
##'   \item{res.all}{A list with all the optimisations results, each element contains
##'                  (almost) the same information as \code{res.best}, e.g.
##'                  \code{res.all[[i]]} contains optimisation results for the i:th
##'                  starting point.}
##'   \item{summary}{A list with parameter estimates and convergence information
##'                  for all starting points.}
##' 
##' @example Rd_examples/Ex_estimate_STmodel.R
##' 
##' @author Johan Lindström
##' @family STmodel methods
##' @family estimateSTmodel methods
##' @method estimate STmodel
##' @export
estimate.STmodel <- function(object, x, x.fixed=NULL, type="p",
                             h=1e-3, diff.type=1, hessian.all=FALSE,
                             lower=-15, upper=15, method="L-BFGS-B",
                             control=list(trace=3, maxit=1000),
                             restart=0, ...){
  ##check class belonging
  stCheckClass(object, "STmodel", name="object")
  
  ##get size of the models
  dimensions <- loglikeSTdim(object)
  
  ##ensure lower case
  type <- tolower(type) 
  ##first check that type is valid
  stCheckType(type)
  
  ##Second check the input x
  x <- as.matrix(x)
  ##if x has length one we need to construct a matrix of initial values
  if( length(x)==1 ){
    x <- matrix(seq(-5,5,length.out=x), dimensions$nparam.cov, x, byrow=TRUE)
  }
  ##check that starting point is valid
  tmp <- stCheckX(x, x.fixed, dimensions, type, object)
  x.all <- tmp$x.all
  x <- tmp$x
  x.fixed <- tmp$x.fixed

  ##Default values for control
  control <- defaultList(control, eval(formals(estimate.STmodel)$control) )
                                    
  ##define local version of gradient function, fixing h and diff.type
  loglikeSTGrad.loc <- function(x, STmodel, type, x.fixed){
    loglikeSTGrad(x, STmodel, type, x.fixed, h=h, diff.type=diff.type)
  }##function loglikeSTGrad.loc

  ##attempt to fit the model for each of the provided starting values.
  res <- as.list( rep(NA, dim(x)[2]) )
  ##vector with convergence information and optimal values
  conv <- rep(FALSE, dim(x)[2])
  value <- rep(NA, dim(x)[2])

  ##make sure that likelihood evaluates
  err <- tryCatch( loglikeST(x[,1], STmodel=object, type=type,
                             x.fixed=x.fixed), silent=TRUE )
  if( inherits(err,"try-error") ){
    stop( paste("log-likelihood fails at first starting point with:\n",
                err[[1]]) )
  }
  
  ##ensure that we are doing maximisation
  control$fnscale <- -1
  ##loop over starting values
  for(i in 1:dim(x)[2]){
    if( control$trace!=0 ){
      message( paste("Optimisation using starting value ",
                     i, "/", dim(x)[2], sep="") )
    }
    x.start <- x[,i]
    i.restart <- 0
    while(i.restart<=restart && !conv[i]){
      try( res[[i]] <- optim(x.start, loglikeST, gr=loglikeSTGrad.loc,
                             STmodel=object, type=type, x.fixed=x.fixed,
                             method=method, control=control, hessian=TRUE,
                             lower=lower, upper=upper), silent=TRUE)
      if( all( !is.na(res[[i]]) ) ){
        ##optim done, let's see if we've converged and update starting point
        x.start <- res[[i]]$par
        ##compute convergence criteria
        conv[i] <- (res[[i]]$convergence==0 &&
                    all(eigen(res[[i]]$hessian)$value < -1e-10))
      }else{
        ##error occured in optim, break
        break
      }
      ##increase counter
      i.restart <- i.restart+1
    }##while(i.restart<=restart && !optim.done)
    if( control$trace!=0 ){
      message("") ##spacing
    }
    
    ##has optimisation finished with out failing?
    if( all( !is.na(res[[i]]) ) ){
      ##extract ML-value
      value[i] <- res[[i]]$value

      ##add convergence and initial parameters
      res[[i]]$conv <- conv[i]
      res[[i]]$par.cov <- data.frame(par=double(dimensions$nparam.cov), sd=NA,
                                     fixed=double(dimensions$nparam.cov),
                                     init=double(dimensions$nparam.cov),
                                     tstat=double(dimensions$nparam.cov))
      res[[i]]$par.all <- data.frame(par=double(dimensions$nparam), sd=NA,
                                     fixed=double(dimensions$nparam),
                                     init=double(dimensions$nparam),
                                     tstat=double(dimensions$nparam))
     
      ##add standard deviations
      suppressWarnings( par.sd <- sqrt(-diag(solve(res[[i]]$hessian))) )
      ##initial value
      if( type!="f" ){
        par.type <- "par.cov"
      }else{
        par.type <- "par.all"
      }
      ##parameters
      res[[i]][[par.type]]$init <- x.all[,i]
      res[[i]][[par.type]]$par <- res[[i]][[par.type]]$fixed <- x.fixed
      res[[i]][[par.type]]$par[is.na(x.fixed)] <- res[[i]]$par
      ##standard error
      res[[i]][[par.type]]$sd[is.na(x.fixed)] <- par.sd
      
      if( type!="f" ){
        ##compute regression parameters
        tmp <- predict(object, res[[i]]$par.cov$par, only.pars=TRUE,
                       pred.var=FALSE, type=type)$pars
        res[[i]]$par.all$par <- c(tmp$gamma.E, tmp$alpha.E, res[[i]]$par.cov$par)
        N.reg <- length(tmp$gamma.E)+length(tmp$alpha.E)
        res[[i]]$par.all$sd <- c(rep(NA,N.reg), res[[i]]$par.cov$sd)
        res[[i]]$par.all$init <- c(rep(NA,N.reg), x.all[,i])
        res[[i]]$par.all$fixed <- c(rep(NA,N.reg), x.fixed)
      }else{
        ##all the covariance parameters
        I <- (dimensions$nparam-dimensions$nparam.cov+1):dimensions$nparam
        res[[i]]$par.cov <- res[[i]]$par.all[I,,drop=FALSE]
      }
      ##t-statistic
      res[[i]]$par.cov$tstat <- res[[i]]$par.cov$par / res[[i]]$par.cov$sd
      res[[i]]$par.all$tstat <- res[[i]]$par.all$par / res[[i]]$par.all$sd
      
      ##add names to the variables.
      rownames(res[[i]]$par.all) <- loglikeSTnames(object, all=TRUE)
      rownames(res[[i]]$par.cov) <- loglikeSTnames(object, all=FALSE)

      if( type!="f" ){
        names(res[[i]]$par) <- loglikeSTnames(object, all=FALSE)[is.na(x.fixed)]
      }else{
        names(res[[i]]$par) <- loglikeSTnames(object, all=TRUE)[is.na(x.fixed)]
      }
    }##if(all(!is.na(res[[i]])))
  }##for(i in 1:dim(x)[2])

  if( all(is.na(res)) ){
    stop("All optimisations failed, consider trying different starting values.")
  }

  ##extract summaries of the optimisations
  status <- data.frame(value=value, convergence=logical(length(res)),
                       conv=(conv==1))
  par.cov <- matrix(NA, dimensions$nparam.cov, length(res))
  par.all <- matrix(NA, dimensions$nparam, length(res))
  for(i in 1:length(res)){
    if( all(!is.na(res[[i]])) ){
      status$convergence[i] <- res[[i]]$convergence==0
      par.cov[,i] <- res[[i]]$par.cov$par
      par.all[,i] <- res[[i]]$par.all$par
    }
  }
  ##add names to the summaries
  rownames(par.all) <- loglikeSTnames(object, all=TRUE)
  rownames(par.cov) <- loglikeSTnames(object, all=FALSE)
  
  ##pick out the converged option with the best value
  Ind.overall <- which.max(value)
  if(any(conv==TRUE)){
    ##mark no-converged values as NA to avoid picking these
    value[!conv] <- NA
  }
  ##extract the best value
  Ind <- which.max(value)
  res.best <- res[[Ind]]
  
  ##collect status results
  summary <- list(status=status, par.all=par.all, par.cov=par.cov, x.fixed=x.fixed)
  
  if(hessian.all==TRUE){
    if(type!="f"){
      x.fixed <- res.best$par.all$fixed
      x <- res.best$par.all$par[ is.na(x.fixed) ]
      res.best$hessian.all <- loglikeSTHessian(x, object, type="f",
                                               x.fixed=x.fixed, h=h)
      ##standard error
      suppressWarnings( par.sd <- sqrt(-diag(solve(res.best$hessian.all))) )
      res.best$par.all$sd <- NA
      res.best$par.all$sd[ is.na(x.fixed) ] <- par.sd
      ##update t-statistic for the best result, all parameters
      res.best$par.all$tstat <- res.best$par.all$par/res.best$par.all$sd
    }else{
      ##replicate hessian for all parameters so output is consistent
      res.best$hessian.all <- res.best$hessian
    }
  }##if(hessian.all==TRUE)
  
  ##return result
  out <- list(res.best=res.best, res.all=res, summary=summary)
  class(out) <- "estimateSTmodel"
  
  return( out )
}##function estimate.STmodel

####################################
## S3 methods for estimateSTmodel ##
####################################
##' \code{\link[base:print]{print}} method for class \code{estimateSTmodel}.
##'
##' @title Print details for \code{estimateSTmodel} object
##' @param x \code{estimateSTmodel} object to print information for.
##' @param ... Ignored additional arguments.
##' @return Nothing
##'
##' @author Johan Lindström
##'
##' @examples
##'   ##load data
##'   data(est.mesa.model)
##'   print(est.mesa.model)
##' 
##' @family estimateSTmodel methods
##' @method print estimateSTmodel
##' @export
print.estimateSTmodel <- function(x, ...){
  ##check class belonging
  stCheckClass(x, "estimateSTmodel", name="x")

  N.opt <- length(x$res.all)
  N.conv <- sum(x$summary$status$conv)
  N.failed <- sum( is.na(x$summary$status$value) )
  I.best <- which.max( x$summary$status$value )
  if( N.conv>0 ){
    tmp <- x$summary$status$value
    tmp[!x$summary$status$conv] <- NA
    I.best.conv <- which.max(tmp)
  }

  cat("Optimisation for STmodel with", N.opt, "starting points.\n")
  cat("  Results:", N.conv, "converged,", N.opt-N.failed-N.conv,
      "not converged,", N.failed,"failed.\n")
  cat("  Best result for starting point", I.best)
  if( x$summary$status$conv[I.best] ){
    cat(", optimisation has converged\n")
  }else{
    cat(", optimisation has NOT converged\n")
    if( N.conv>0 ){
      cat("  Best converged result for starting point", I.best.conv, "\n")
    }
  }
  cat("\n")

  ##fixed points
  if( all(is.na(x$summary$x.fixed)) ){
    cat("No fixed parameters.\n")
  }else{
    cat("Fixed parameters:\n")
    print(x$summary$x.fixed[!is.na(x$summary$x.fixed)])
  }
  
  ##estimated parameters
  cat("\nEstimated parameters for all starting point(s):\n")
  print(x$summary$par.all)
  cat("\nFunction value(s):\n")
  print(x$summary$status$value)
  cat("\n")

  return(invisible())
}##function print.estimateSTmodel

##' \code{\link[stats:coef]{coef}} method for class \code{estimateSTmodel}.
##'
##' @title Returns estimated parameters (and uncertaintes)
##' @param object \code{estimateSTmodel} object from which to extract estimated
##'   parameters.
##' @param pars One of "cov", "reg", "all"; which parameters to extract.
##' @param ... Ignored additional arguments.
##' @return Estimated parameters.
##'
##' @author Johan Lindström
##'
##' @examples
##'   ##load data
##'   data(est.mesa.model)
##'   ##extract all parameters
##'   coef(est.mesa.model)
##'   ##extract only covariance parameters
##'   coef(est.mesa.model, pars="cov")
##' 
##' @family estimateSTmodel methods
##' @importFrom stats coef
##' @method coef estimateSTmodel
##' @export
coef.estimateSTmodel <- function(object, pars=c("all", "cov", "reg"), ...){
  ##check class belonging
  stCheckClass(object, "estimateSTmodel", name="object")

  pars <- match.arg(pars)

  ##pick which parameters
  if( pars=="cov" ){
    res <- object$res.best$par.cov
  }else if( pars=="reg" ){
    res <- object$res.best$par.all[1:(dim(object$res.best$par.all)[1] -
                                       dim(object$res.best$par.cov)[1]),,drop=FALSE]
  }else if( pars=="all" ){
    res <- object$res.best$par.all
  }
  return(res)
}##function coef.estimateSTmodel

#############################
## S3 methods for estimate ##
#############################

##' @rdname estimate.STmodel
##' @export
estimate <- function(object, x, ...){
  UseMethod("estimate")
}
