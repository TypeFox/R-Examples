###################################
## Functions for crossvalidation ##
###################################
##Functions in this file:
## estimateCV.STmodel          EX:ok
## estimateCV                  EX:S3 method
## print.estCVSTmodel          EX:ok
## summary.estCVSTmodel        EX:ok
## print.summary.estCVSTmodel  EX:with summary.estCVSTmodel
## coef.estCVSTmodel           EX:ok
## boxplot.estCVSTmodel        EX:ok

##' Functions that perform cross-validated parameter estimation and prediction
##' for the spatio-temporal model.
##' 
##' For \code{predictCV.STmodel} the parameters used to compute predictions for the left
##' out observations can be either a single vector or a matrix.
##' For a single vector the same parameter values will be used for all
##' cross-validation predictions; for a matrix the parameters in \code{x[,i]}
##' will be used for the predictions of the i:th cross-validation set (i.e. for
##' \code{Ind.cv[,i]}). Suitable matrices are provided in the output from
##' \code{estimateCV.STmodel}.
##' 
##' The cross-validation groups are given by \code{Ind.cv}. \code{Ind.cv} should
##' be either a (number of observations) - by - (groups) logical matrix or an
##' \emph{integer valued} vector with length equal to (number of observations).
##' If a matrix then each column defines a cross-validation set with the
##' \code{TRUE} values marking the observations to be left out. If a vector then
##' \code{1}:s denote observations to be dropped in the first cross-validation
##' set, \code{2}:s observations to be dropped in the second set, etc.
##' Observations marked by values \code{<=0} are never dropped. See
##' \code{\link{createCV}} for details.
##'
##' @title Cross-Validated Estimation and Prediction
##' 
##' @param object \code{STmodel} object for which to perform cross-validation.
##' @param x Either a vector or matrix of starting point(s) for the optimisation,
##'   see \code{\link{estimate.STmodel}}; or a matrix with parameters, the i:th
##'   row being used for prediction of the i:th cross-validation set. For
##'   prediction either a \code{estCVSTmodel} or \code{estimateSTmodel} object,
##'   results from \code{\link{estimateCV.STmodel}} or
##'   \code{\link{estimate.STmodel}}, can be used. 
##' @param Ind.cv \code{Ind.cv} defines the cross-validation scheme.  Either a
##'   (number or observations) - by - (groups) logical matrix or an \emph{integer
##'   valued} vector with length equal to (number or observations). For
##'   \code{predictCV.STmodel} \code{Ind.cv} can be infered from \code{x} if
##'   \code{x} is a \code{estCVSTmodel} object
##'   See further \code{\link{createCV}}.
##' @param control A list of control parameters for the optimisation.
##'   See \code{\link[stats:optim]{optim}} for details; setting \code{trace}=0
##'   eliminates all ouput.
##' @param verbose.res A \code{TRUE}/\code{FALSE} variable indicating if full
##'   results from \code{\link{estimate.STmodel}} for each CV group should be
##'   returned; defaults to \code{FALSE}
##' @param ... All additional parameters for \code{\link{estimate.STmodel}}
##'   or \code{\link{predict.STmodel}}.
##'   For \code{\link{predict.STmodel}} a number of parameters are set in
##'   \code{predictCV.STmodel} and can \strong{NOT} be overriden, these are
##'   \code{nugget.unobs}, \code{only.pars=FALSE}, and
##'   \code{combine.data=FALSE}. 
##' 
##' @return Either a \code{estCVSTmodel} object with elements:
##'   \item{status}{Data.frame with convergence information and best function
##'                 value for each cross-validation group.}
##'   \item{Ind.cv}{The cross-validation grouping.}
##'   \item{x.fixed}{Fixed parameters in the estimation, see
##'                  \code{\link{estimate.STmodel}}.}
##'   \item{x.init}{Matrix of inital values used, i.e. \code{x} from the input.}
##'   \item{par.all, par.cov}{Matrices with estimated parameters for each
##'                           cross-validation group.}
##'   \item{par.all.sd, par.cov.sd}{Standard deviations computed from the
##'     Hessian/information matrix for set of estimated parameters.}
##'   \item{res.all}{Estimation results for each cross-validation group,
##'                  contains the output from the \code{\link{estimate.STmodel}}
##'                  calls, only included if \code{verbose.res=TRUE}.}
##' Or a \code{predCVSTmodel} object with elements:
##'   \item{opts}{Copy of the \code{opts} field in the output from
##'               \code{\link{predict.STmodel}}.}
##'   \item{Ind.cv}{The cross-validation grouping.}
##'   \item{pred.obs}{A data.frame with a copy of observations from
##'                   \code{object$obs}, predictions (for different model
##'                   components), variances, and residuals. Variance field will
##'                   be missing if \code{pred.var=FALSE}.}
##'   \item{pred.all}{A list with time-by-location data.frames containing
##'                   predictions and variances for all space-time locations as
##'                   well as predictions and variances for the
##'                   beta-fields. Unobserved points are \code{NA} for the
##'                   option \code{only.obs=TRUE}.}
##' 
##' @example Rd_examples/Ex_estimateCV_STmodel.R
##'
##' @author Johan Lindström
##' @family STmodel methods
##' @family cross-validation functions
##' @family estCVSTmodel methods
##' @method estimateCV STmodel
##' @export
estimateCV.STmodel <- function(object, x, Ind.cv, control=list(trace=3),
                               verbose.res=FALSE, ...){
  ##check class belonging
  stCheckClass(object, "STmodel", name="object")
  ##check cross-validation groups
  Ind.cv <- stCheckInternalCV(Ind.cv)
  
  ##Default values for control
  control <- defaultList(control, eval(formals(estimate.STmodel)$control) )
  
  ##ensure that Ind.cv is a matrix
  Ind.cv <- as.matrix(Ind.cv)
  if( dim(Ind.cv)[2]==1 ){
    N.CV.sets <- max(Ind.cv,na.rm=TRUE)
  }else{
    N.CV.sets <- dim(Ind.cv)[2]
  }
  
  ##get size of the models
  dimensions <- loglikeSTdim(object)
  
  res <- list()
  for(i in 1:N.CV.sets){
    if( dim(Ind.cv)[2]==1 ){
      Ind.current <- Ind.cv==i
    }else{
      Ind.current <- as.logical( Ind.cv[,i] )
    }
    ##create data matrices that contains observations
    object.aux <- dropObservations(object, Ind.current)
    
    ##lets estimate parameters for this set
    if( control$trace!=0 ){
      message( "\n***************************")
      message( paste("Estimation of CV-set ", i, "/", N.CV.sets, sep="") )
    }
    res[[i]] <- estimate(object=object.aux, x=x, control=control, ...)
  }##for(i in 1:dim(Ind.cv)[2])

  ##status of optimisations
  status <- data.frame(value=sapply(res, function(x){x$res.best$value}),
                       convergence=sapply(res, function(x){x$res.best$convergence==0}),
                       conv=sapply(res, function(x){x$res.best$conv}))
  ##add hessian eigen-values to status
  tmp <- sapply(res, function(x){ range( -eigen(x$res.best$hessian)$value) })
  status$eigen.min <- tmp[1,]
  if( is.null( res[[1]]$res.best$hessian.all) ){
    status$eigen.all.min <- NA
  }else{
    tmp <- sapply(res,
                  function(x){ range( -eigen(x$res.best$hessian.all)$value) })
    status$eigen.all.min <- tmp[1,]
  }
  
  ##matrices with the estimates parameters (accounting for x.fixed)
  par.cov <- matrix(NA, dim(res[[1]]$res.best$par.cov)[1], N.CV.sets)
  par.cov.sd <- matrix(NA, dim(res[[1]]$res.best$par.cov)[1], N.CV.sets)
  par.all <- matrix(NA, dim(res[[1]]$res.best$par.all)[1], N.CV.sets)
  par.all.sd <- matrix(NA, dim(res[[1]]$res.best$par.all)[1], N.CV.sets)

  ##best estimates for each CV-group
  for(i in 1:N.CV.sets){
    par.cov[,i] <- res[[i]]$res.best$par.cov$par
    par.cov.sd[,i] <- res[[i]]$res.best$par.cov$sd
    par.all[,i] <- res[[i]]$res.best$par.all$par
    par.all.sd[,i] <- res[[i]]$res.best$par.all$sd
  }
  ##all initial values
  x.init <- sapply(res[[1]]$res.all, function(x){ x$par.all$init})
  ##also return a list with all optimisation results?
  if( verbose.res ){
    res.all <- res
  }else{
    res.all <- NULL
  }
  
  rownames(par.cov) <- rownames(res[[1]]$res.best$par.cov)
  rownames(par.cov.sd) <- rownames(res[[1]]$res.best$par.cov)
  rownames(par.all) <- rownames(res[[1]]$res.best$par.all)
  rownames(par.all.sd) <- rownames(res[[1]]$res.best$par.all)
  rownames(x.init) <- rownames(res[[1]]$res.best$par.all)
  ##drop NA's from x.init (occurs if only covariance parameters where given)
  x.init <- x.init[!apply(is.na(x.init),1,all),,drop=FALSE]
  ##Return the estimated parameters from the cross-validation
  out <- list(par.cov=par.cov, par.cov.sd=par.cov.sd,
              par.all=par.all, par.all.sd=par.all.sd,
              res.all=res.all, status=status, Ind.cv=Ind.cv,
              x.fixed=res[[1]]$summary$x.fixed,
              x.init=x.init)
  class(out) <- "estCVSTmodel"

  return( out )
}##function estimateCV.STmodel

#######################################
## General S3 methods for estimateCV ##
#######################################

##' @rdname estimateCV.STmodel
##' @export
estimateCV <- function(object, x, Ind.cv, ...){
  UseMethod("estimateCV")
}

#################################
## S3-METHODS FOR estCVSTmodel ##
#################################

##' \code{\link[base:print]{print}} method for class \code{estCVSTmodel}.
##'
##' @title Print details for \code{estCVSTmodel} object
##' @param x \code{estCVSTmodel} object to print information for.
##' @param ... Ignored additional arguments.
##' @return Nothing
##' 
##' @examples
##' ##load some data
##' data(est.cv.mesa)
##' ##print basic information for the CV-predictions
##' print(est.cv.mesa)
##' 
##' @author Johan Lindström
##' 
##' @family estCVSTmodel methods
##' @method print estCVSTmodel
##' @export
print.estCVSTmodel <- function(x, ...){
  ##check class belonging
  stCheckClass(x, "estCVSTmodel", name="x")

  N.cv <- dim(x$status)[1]
  N.opt <- dim(x$x.init)[2]
  N.conv <- sum(x$status$conv)

  cat("Cross-validation parameter estimation for STmodel\n")
  cat("  with", N.cv, "CV-groups and", N.opt, "starting points.\n")
  cat("  Results:", N.conv, "converged,", N.cv-N.conv, "not converged.\n")
  cat("\n")

  ##fixed points
  if( all(is.na(x$x.fixed)) ){
    cat("No fixed parameters.\n")
  }else{
    cat("Fixed parameters:\n")
    print(x$x.fixed[!is.na(x$x.fixed)])
  }
  cat("\n")

  cat("Estimated function values and convergence info:\n")
  print(x$status)
  cat("\n")

  return(invisible())
}##function print.estCVSTmodel

##' \code{\link[base:summary]{summary}} method for class \code{estCVSTmodel}.
##'
##' @title Computes summary details for \code{estCVSTmodel} object
##' @param object \code{estCVSTmodel} object to compute summary information for.
##' @param ... Ignored additional arguments. 
##' @return A \code{summary.estCVSTmodel} object.
##' 
##' @examples
##' ##load some data
##' data(est.cv.mesa)
##' ##print basic information for the CV-predictions
##' summary(est.cv.mesa)
##'
##' @author Johan Lindström
##' 
##' @family estCVSTmodel methods
##' @method summary estCVSTmodel
##' @export
summary.estCVSTmodel <- function(object, ...){
  ##check class belonging
  stCheckClass(object, "estCVSTmodel", name="object")
  
  ##allocate output object
  out <- list()
  ##compute various summaries, estimated parameters
  out$value <- summary(object$status$value)
  out$par.cov <- summary( t(object$par.cov) )
  out$par.all <- summary( t(object$par.all) )

  ##return the object
  class(out) <- "summary.estCVSTmodel"
  return(out)
}##function summary.estCVSTmodel

##' \code{\link[base:print]{print}} method for class \code{summary.estCVSTmodel}.
##'
##' @title Print details for \code{summary.estCVSTmodel} object
##' @param x \code{summary.estCVSTmodel} object to print information for.
##' @param ... Additional arguments, passed to
##'   \code{\link[base:print]{print.table}}.
##' @return Nothing
##'
##' @author Johan Lindström
##' 
##' @family estCVSTmodel methods
##' @method print summary.estCVSTmodel
##' @export
print.summary.estCVSTmodel <- function(x, ...){
  ##check class belonging
  stCheckClass(x, "summary.estCVSTmodel", name="x")

  ##print data
  ##estimated parameters
  cat("Cross-validation parameter estimation for STmodel:\n")
  cat("Summary for function values:\n")
  print( x$value, ...)
  cat("\n")

  cat("Summary for all parameters:\n")
  print( x$par.all, ...)
  cat("\n")
  
  return(invisible())
}##function print.summary.estCVSTmodel

##' \code{\link[stats:coef]{coef}} method for class \code{estCVSTmodel}.
##'
##' @title Returns estimated parameters for each CV-group.
##' @param object \code{estCVSTmodel} object from which to extract estimated
##'   parameters.
##' @param pars One of "cov", "reg", "all"; which parameters to extract.
##' @param ... Ignored additional arguments.
##' @return Nothing
##'
##' @author Johan Lindström
##'
##' @examples
##'   ##load data
##'   data(est.cv.mesa)
##'   ##extract all parameters
##'   coef(est.cv.mesa)
##'   ##extract only covariance parameters
##'   coef(est.cv.mesa, pars="cov")
##' 
##' @family estCVSTmodel methods
##' @importFrom stats coef
##' @method coef estCVSTmodel
##' @export
coef.estCVSTmodel <- function(object, pars=c("all", "cov", "reg"), ...){
  ##check class belonging
  stCheckClass(object, "estCVSTmodel", name="object")

  pars <- match.arg(pars)

  ##pick which parameters
  if( pars=="cov" ){
    res <- object$par.cov
  }else if( pars=="reg" ){
    res <- object$par.all[1:(dim(object$par.all)[1] -
                             dim(object$par.cov)[1]),,drop=FALSE]
  }else if( pars=="all" ){
    res <- object$par.all
  }
  return(res)
}##function coef.estimateSTmodel

##' \code{\link[graphics:boxplot]{boxplot}} method for class \code{estCVSTmodel}.
##' 
##' @title Boxplots for \code{estCVSTmodel} object
##' @param x \code{estCVSTmodel}/\code{STmodel} object to boxplot.
##' @param plot.type One of "cov", "reg", "all"; should we boxplot covariance,
##'   regression or all parameter estimates.
##' @param ... Additional parameters passed to
##'   \code{\link[graphics:boxplot]{boxplot}}.
##' @return Nothing
##'
##' @example Rd_examples/Ex_boxplot_estCVSTmodel.R
##'
##' @author Johan Lindström
##' 
##' @family estCVSTmodel methods
##' @importFrom graphics boxplot
##' @method boxplot estCVSTmodel
##' @export
boxplot.estCVSTmodel <- function(x, plot.type=c("cov", "reg", "all"), ...){
  ##check class belonging
  stCheckClass(x, "estCVSTmodel", name="x")

  plot.type <- match.arg(plot.type)

  ##pick which parameters
  if( plot.type=="cov" ){
    pars <- x$par.cov
  }else if( plot.type=="reg" ){
    pars <- x$par.all[1:(dim(x$par.all)[1]-dim(x$par.cov)[1]),,drop=FALSE]
  }else if( plot.type=="all" ){
    pars <- x$par.all
  }

  ##boxplot for parameters
  boxplot(t(pars), ...)
  
  return(invisible())
}##function boxplot.estCVSTmodel
