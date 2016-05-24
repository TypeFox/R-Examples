########################
## Functions for MCMC ##
########################
##Functions in this file:
## MCMC.STmodel               EX:OK
## print.mcmcSTmodel          EX:OK
## summary.mcmcSTmodel        EX:OK
## print.summary.mcmcSTmodel  EX:with summary.estCVSTmodel
## plot.mcmcSTmodel           EX:OK
## density.mcmcSTmodel        EX:OK
## plot.density.mcmcSTmodel   EX:with density.mcmcSTmodel

##' Estimates parameters and parameter uncertainties for the spatio-temporal
##' model using a Metropolis-Hastings based Markov Chain Monte Carlo (MCMC)
##' algorithm.
##' \cr
##' The function runs uses a Metropolis-Hastings algorithm (Hastings, 1970) to
##' sample from the parameters of the spatio-temporal model, assuming flat
##' priors for all the parameters (flat on the log-scale for the covariance
##' parameters).
##'
##' At each iteration of the MCMC new parameters are proposed using a
##' random-walk with a proposal covariance matrix. The proposal matrix is
##' determined as:
##' \describe{
##'   \item{1}{If \code{Sigma.prop} is given then this is used.}
##'   \item{2}{If \code{Sigma.prop=NULL} then we follow Roberts et.al. (1997) and
##'            compute \cr
##'            \code{c <- 2.38*2.38/dim(Hessian.prop)[1]}\cr 
##'            \code{Sigma.prop <- -c*solve(Hessian.prop)}.}
##'   \item{3}{If both \code{Sigma.prop=NULL} and \code{Hessian.prop=NULL} then the
##'            Hessian is computed using \code{\link{loglikeSTHessian}} and
##'            \code{Sigma.prop} is computed according to point 2.}
##' }
##' The resulting proposal matrix is checked to ensure that it is positive
##' definite before proceeding, \cr
##' \code{all(eigen(Sigma.prop)$value > 1e-10)}.
##'
##' @title MCMC Inference of Parameters in the Spatio-Temporal Model
##' 
##' @param object \code{STmodel} for which to run MCMC.
##' @param x Point at which to start the MCMC. Could be either only
##'   \emph{log}-covariance parameters or regression and \emph{log}-covariance
##'   parameters. If regression parameters are given but not needed they are
##'   dropped, if they are needed but not given they are inferred by calling \cr
##'   \code{\link{predict.STmodel}} with \code{only.pars=TRUE}.
##' @param x.fixed Vector with parameter to be held fixed; parameters marked as
##'   \code{NA} will still be estimated.
##' @param type A single character indicating the type of log-likelihood to
##'   compute. Valid options are "f" or "r", for \emph{full}, or \emph{restricted 
##'   maximum likelihood} (REML). Since \emph{profile} is not a proper
##'   likelihood \code{type="p"} will revert (with a warning) to using the
##'   \emph{full} log-likelihood.
##' @param N Number of MCMC iterations to run.
##' @param Hessian.prop Hessian (information) matrix for the log-likelihood, can
##'   be used to create a proposal matrix for the MCMC.
##' @param Sigma.prop Proposal matrix for the MCMC.
##' @param info Outputs status information every info:th iteration.
##'   If \code{info=0} no output.
##' @param ... ignored additional arguments.
##' 
##' @return \code{mcmcSTmodel} object with elements:
##'   \item{par}{A \code{N} - by - (number of parameters) matrix with
##'              trajectories of the parameters.}
##'   \item{log.like}{A vector of length \code{N} with the log-likelihood values
##'                    at each iteration.}
##'   \item{acceptance}{A vector of length \code{N} with the acceptance
##'                     probability for each iteration.}
##'   \item{Sigma.prop, chol.prop}{Proposal matrix and it's Choleskey factor.}
##'   \item{x.fixed}{Any fixed parameters.}
##' 
##' @example Rd_examples/Ex_MCMC_mesa_model.R
##' 
##' @author Johan Lindström
##' @family STmodel methods
##' @family mcmcSTmodel methods
##' @method MCMC STmodel
##' @export
MCMC.STmodel <- function(object, x, x.fixed=NULL, type="f", N=1000,
                         Hessian.prop=NULL, Sigma.prop=NULL,
                         info=min(ceiling(N/50),100),... ){
  ##check class belonging
  stCheckClass(object, "STmodel", name="object")
  
  ##get size of the models
  dimensions <- loglikeSTdim(object)
  
  ##ensure lower case
  type <- tolower(type) 
  ##first check that type is valid
  stCheckType(type)
  if( type=="p" ){
    warning("Profile likelihood not proper, using the full parameter set (type=\"f\")",
            immediate.=TRUE)
    type <- "f"
  }
  ##check that starting point is valid
  x <- matrix(x, ncol=1)
  tmp <- stCheckX(x, x.fixed, dimensions, type, object)
  x.all <- tmp$x.all
  x <- tmp$x
  x.fixed <- tmp$x.fixed
  
  if( is.null(Sigma.prop) || any(is.na(Sigma.prop)) ){ ##no proposal given
    ##see if we need to compute the Hessian, either due to missing proposal
    if( is.null(Hessian.prop) || any(is.na(Hessian.prop)) ){
      if( info!=0 ){
        message("Computing Hessian")
      }
      Hessian.prop <- loglikeSTHessian(x, object, type=type, x.fixed=x.fixed)
    }else if( any(dim(Hessian.prop)!=length(x)) ){
      ##or due to size missmatch
      warning("Size missmatch, recomputing Hessian", immediate.=TRUE)
      Hessian.prop <- loglikeSTHessian(x, object, type=type, x.fixed=x.fixed)
    }
    ##check if Hessian is PSD
    if( !all(eigen(Hessian.prop)$value < -1e-10) ){
      stop("Hessian not negative definite (some eigs > -1e-10)")
    }
    ##compute proposal matrix from the Hessian ( inv(H)*2.38/d )
    Sigma.prop <- -solve(Hessian.prop)*2.38*2.38/dim(Hessian.prop)[1]
  }
  
  ##check size of given proposal
  if( any(dim(Sigma.prop)!=length(x)) ){
    stop( paste("'Sigma.prop' does not match parameter vector, should be ",
                length(x), "-by-", length(x), sep="") )
  }
  ##check PSD
  if( !all(eigen(Sigma.prop)$value > 1e-10) ){
    stop("Proposal matrix not positive definite (some eigs < 1e-10)")
  }
  ##compute choleskey
  S <- t( chol(Sigma.prop) )

  ##allocate storage
  l <- double(N)
  alpha <- double(N)
  par <- matrix(NA, N, length(x.fixed))
  ##parameter names
  colnames(par) <- names(x.fixed)
  colnames(S) <- rownames(S) <- colnames(par)[is.na(x.fixed)]

  ##initial values
  par[1,] <- x.fixed
  par[1,is.na(x.fixed)] <- x
  l[1] <- loglikeST(par[1,], object, type=type)
  par.new <- par[1,]
  ##run MCMC
  for(i in 2:N){
    ##propose new value
    par.new[is.na(x.fixed)] <- par[i-1,is.na(x.fixed)] + S %*% rnorm(dim(S)[1])
    ##evaluate likelihood for the new value
    l[i] <- loglikeST(par.new, object, type=type)
    alpha[i] <- min(exp(l[i]-l[i-1]),1)
    if( runif(1)<alpha[i] ){ #accept
      par[i,] <- par.new
    }else{#reject
      par[i,] <- par[i-1,]
      l[i] <- l[i-1]
    }
    if(info!=0 && (i %% info)==0)
      message( sprintf("MCMC %d/%d",i,N) )
  }#for(i in 2:N)
  
  ##return the results
  out <- list(par=par, log.like=l, acceptance=alpha, Sigma.prop=Sigma.prop,
              chol.prop=S, x.fixed=x.fixed)
  class(out) <- "mcmcSTmodel"
  return( out )
}#MCMC.STmodel

#################################
## General S3 methods for MCMC ##
#################################
##' @rdname MCMC.STmodel
##' @export
MCMC <- function(object, ...){
  UseMethod("MCMC")
}

################################
## S3 methods for mcmcSTmodel ##
################################
##' \code{\link[base:print]{print}} method for class \code{mcmcSTmodel}.
##'
##' @title Print details for \code{mcmcSTmodel} object
##' @param x \code{mcmcSTmodel} object to print information for.
##' @param ... Ignored additional arguments.
##' @return Nothing
##'
##' @examples
##'   ##load data
##'   data(MCMC.mesa.model)
##'   print(MCMC.mesa.model)
##' 
##' @author Johan Lindström
##'
##' @family mcmcSTmodel methods
##' @method print mcmcSTmodel
##' @export
print.mcmcSTmodel <- function(x, ...){
  ##check class belonging
  stCheckClass(x, "mcmcSTmodel", name="x")

  N <- length(x$l)
  N.pars <- dim(x$par)[2]
  p <- mean(x$acceptance)
  
  cat("MCMC for STmodel, results over", N, "iterations.\n")
  cat(" ", N.pars, "unknown parameters\n")
  cat("  Acceptance rate:", p, "\n")
  cat("\n")

  ##fixed points
  if( all(is.na(x$x.fixed)) ){
    cat("No fixed parameters.\n")
  }else{
    cat("Fixed parameters:\n")
    print(x$x.fixed[!is.na(x$x.fixed)])
  }
  cat("\n")

  return(invisible())
}##function print.estimateSTmodel

##' \code{\link[base:summary]{summary}} method for class \code{mcmcSTmodel}.
##'
##' @title Computes summary details for \code{mcmcSTmodel} object
##' @param object \code{mcmcSTmodel} object to compute summary information for.
##' @param burnIn Number of initial iterations to drop.
##' @param ... Ignored additional arguments. 
##' @return A \code{summary.mcmcSTmodel} object.
##'
##' @examples
##'   ##load data
##'   data(MCMC.mesa.model)
##'   summary(MCMC.mesa.model)
##' 
##' @author Johan Lindström
##' 
##' @family mcmcSTmodel methods
##' @method summary mcmcSTmodel
##' @export
summary.mcmcSTmodel <- function(object, burnIn=0, ...){
  ##check class belonging
  stCheckClass(object, "mcmcSTmodel", name="object")

  ##determine burnIn
  I <- rep(TRUE, length(object$l))
  if( burnIn>=length(object$l) ){
    stop( paste("burnIn larger than number of iterations:", length(object$l)) )
  }
  if( burnIn>0 ){
    I[1:burnIn] <- FALSE
  }
  
  ##allocate output object
  out <- list()
  ##compute various summaries, estimated parameters
  out$N <- length(object$l)
  out$burnIn <- burnIn
  out$l <- summary( object$l[I] )
  out$par <- summary( object$par[I,] )
  out$alpha <- summary( object$acceptance[I] )

  ##return the object
  class(out) <- "summary.mcmcSTmodel"
  return(out)
}##function summary.mcmcSTmodel

##' \code{\link[base:print]{print}} method for class \code{summary.mcmcSTmodel}.
##'
##' @title Print details for \code{summary.mcmcSTmodel} object
##' @param x \code{summary.mcmcSTmodel} object to print information for.
##' @param ... Additional arguments, passed to
##'   \code{\link[base:print]{print.table}}.
##' @return Nothing
##'
##' @author Johan Lindström
##' 
##' @family mcmcSTmodel methods
##' @method print summary.mcmcSTmodel
##' @export
print.summary.mcmcSTmodel <- function(x, ...){
  ##check class belonging
  stCheckClass(x, "summary.mcmcSTmodel", name="x")

  ##print data
  ##estimated parameters
  cat("MCMC for STmodel, results over", x$N-x$burnIn, "iterations.\n")
  cat(x$N, "iterations in total.\n")
  cat(x$burnIn, "iterations discarded as burn in.\n\n")
  
  cat("Summary for log-ikelihood values:\n")
  print( x$l, ...)
  cat("\n")
  cat("Summary for acceptance rate:\n")
  print( x$alpha, ...)
  cat("\n")
  cat("Summary for parameters:\n")
  print( x$par, ...)
  cat("\n")
  
  return(invisible())
}##function print.summary.mcmcSTmodel


##' \code{\link[graphics:plot]{plot}} method for class \code{mcmcSTmodel}.
##'
##' Plots results from \code{\link{MCMC.STmodel}}. Either parameter paths or the
##' log-likelihood for the mcmc simulations. 
##' 
##' @title Plots for an \code{mcmcSTmodel} object
##' @param x \code{mcmcSTmodel} object to plot.
##' @param y Type of plot, options are \code{"like"}, \code{"alpha"}, or
##'   name/index number of a parameter.
##' @param add Add to existing plot using \code{\link[graphics:lines]{lines}}
##' @param main Parameter passed as \code{main} to
##'   \code{\link[graphics:plot]{plot}}, defaults to the parameter-name if not
##'   given. 
##' @param ... Additional parameters passed to \code{\link[graphics:plot]{plot}}
##'   or \code{\link[graphics:lines]{lines}}
##' @return Nothing
##'
##' @example Rd_examples/Ex_plot_mcmcSTmodel.R
##' 
##' @author Johan Lindström
##' 
##' @family mcmcSTmodel methods
##' @method plot mcmcSTmodel
##' @export
plot.mcmcSTmodel <- function(x, y="like", add=FALSE, main=NULL, ...){
  ##check class belonging
  stCheckClass(x, "mcmcSTmodel", name="x")

  ##we have to use y, cast to resonable name
  plot.type <- y

  plot.type.OK <- FALSE ##Assume the worst.
  if( is.numeric(plot.type) ){
    plot.type.OK <- 1<=plot.type && plot.type<=dim(x$par)[2]
  }else if( is.character(plot.type) ){
    plot.type.OK <- plot.type %in% c("like","alpha",colnames(x$par))
  }
  if( !plot.type.OK ){
    stop("Unknown option for 'y'")
  }
  if( is.character(plot.type) && (plot.type %in% colnames(x$par)) ){
    plot.type <- match(plot.type,colnames(x$par))
  }

  if( is.character(plot.type) && plot.type=="like" ){
    y <- x$l
    name <- "log-likelihood"
  }else if( is.character(plot.type) && plot.type=="like" ){
    y <- x$acceptance
    name <- "acceptance rate"
  }else{
    y <- x$par[,plot.type]
    name <- colnames(x$par)[plot.type]
  }
  x <- seq(1,length(y))

  if( missing(main) ){ main <- name }
  
  if( add==FALSE ){
    plot(x, y, main=main, ...)
  }else{
    lines(x, y, ...)
  }
  
  return(invisible())
}##function plot.mcmcSTmodel


##' \code{\link[stats:density]{density}} method for class \code{mcmcSTmodel}.
##'
##' Computes kernel density estimates for the MCMC-parameters; as well as
##' approximate Gaussian densities based on the Fischer-information.
##' 
##' @title Kernel Density Estimation for an \code{mcmcSTmodel} Object
##' @param x \code{mcmcSTmodel} object
##' @param BurnIn Number of initial points to ignore.
##' @param estSTmodel Either a \code{estimateSTmodel} object from
##'   \code{\link{estimate.STmodel}} or a matrix with parameter-estimates and
##'   standard deviations, such as the output from
##'   \code{\link{coef.estimateSTmodel}}. If given as a matrix, it should have
##'   columns named "par" and "sd", and rows named after the parameters.
##' @param ... Additional parameters passed to
##'   \code{\link[stats:density]{density}}.
##' @return List containing density estimate and Gaussian densities for all
##'   model parameters.
##'
##' @example Rd_examples/Ex_density_mcmcSTmodel.R
##' 
##' @author Johan Lindström
##' 
##' @family mcmcSTmodel methods
##' @importFrom stats density
##' @method density mcmcSTmodel
##' @export
density.mcmcSTmodel <- function(x, BurnIn=0, estSTmodel=NULL, ...){
  ##check class belonging
  stCheckClass(x, "mcmcSTmodel", name="x")

  ##extract vector of parameters
  if( !is.null(estSTmodel) ){
    if( inherits(estSTmodel, "estimateSTmodel") ){
      estSTmodel <- coef(estSTmodel)
    }
    if( !all(c("par","sd") %in% colnames(estSTmodel)) ){
      stop("estSTmodel must contain columns named par and sd")
    }
    estSTmodel <- estSTmodel[,c("par","sd")]
    if( !all(rownames(estSTmodel) %in% colnames(x$par)) ){
      stop("estSTmodel must contain rows named after the parameters, i.e. colnames(x$par)")
    }
    estSTmodel <- estSTmodel[colnames(x$par),]
  }else{
    estSTmodel <- matrix(NA,dim(x$par)[2],2)
  }

  ##density estimation
  N <- dim(x$par)[1]
  est.dens <- apply(x$par[(BurnIn+1):N,], 2, density, ...)

  ##Gauss approxmation
  colnames(estSTmodel) <- c("mean", "sd")
  rownames(estSTmodel) <- NULL
  
  ##vector containing results
  res <- vector("list", dim(x$par)[2])
  names(res) <- colnames(x$par)
  ##collect results
  for(i in 1:length(res)){
    res[[i]] <- list(density=est.dens[[i]],
                     approx.gauss = list(mean = estSTmodel[i,"mean"],
                       sd = estSTmodel[i,"sd"]))
  }
  class(res) <- "density.mcmcSTmodel"
  return( res )
}##function density.mcmcSTmodel


##' \code{\link[graphics:plot]{plot}} method for class \code{density.mcmcSTmodel}.
##' Plots results from \code{\link{density.mcmcSTmodel}}.
##' 
##' @title Plots for an \code{density.mcmcSTmodel} object
##' @param x \code{density.mcmcSTmodel} object to plot.
##' @param y Name/index of parameter for which to plot the density.
##' @param add Add to existing plot using \code{\link[graphics:lines]{lines}}.
##' @param norm.col Add the Gaussian density using a line with colour
##'   \code{norm.col}, if \code{norm.col=0} do \emph{not} add the Gaussian.
##' @param main Parameter passed as \code{main} to
##'   \code{\link[stats:plot.density]{plot.density}}, defaults to the
##'   parameter-name if not given.
##' @param ylim Additional parameters passed to
##'   \code{\link[stats:plot.density]{plot.density}}.
##' @param ... Additional parameters passed to
##'   \code{\link[stats:plot.density]{plot.density}} or
##'   \code{\link[graphics:lines]{lines}}.
##' @return Nothing
##'
##' @example Rd_examples/Ex_density_mcmcSTmodel.R
##' 
##' @author Johan Lindström
##' 
##' @family mcmcSTmodel methods
##' @method plot density.mcmcSTmodel
##' @export
plot.density.mcmcSTmodel <- function(x, y=1, add=FALSE, norm.col=0,
                                     main=NULL, ylim=NULL, ...){
  ##check class belonging
  stCheckClass(x, "density.mcmcSTmodel", name="x")

  ##we have to use y, cast to resonable name
  plot.type <- y

  plot.type.OK <- FALSE ##Assume the worst.
  if( is.numeric(plot.type) ){
    plot.type.OK <- 1<=plot.type && plot.type<=length(x)
  }else if( is.character(plot.type) ){
    plot.type.OK <- plot.type %in% names(x)
  }
  if( !plot.type.OK ){
    stop("Unknown option for 'y'")
  }

  if( is.numeric(plot.type) ){
    plot.type <- names(x)[plot.type]
  }
  ##extract density and gauss approx to plot
  dens <- x[[plot.type]]$density
  approx <- x[[plot.type]]$approx.gauss
  name <- plot.type

  if( norm.col!=0 ){
    ##Gaussian approximation
    xd <- dens$x
    yd <- dnorm(xd, mean=approx$mean, sd=approx$sd)
  }else{
    yd <- NULL
  }

  if( missing(main) ){ main <- name }

  if( is.null(ylim) ){
    ylim <- range(c(dens$y,yd))
  }else{
    ylim <- range(ylim)
  }

  ##plot density kernel estimate
  if( add==FALSE ){
    plot(dens, main=main, ylim=ylim, ...)
  }else{
    lines(dens$x, dens$y, ...)
  }
  if( norm.col!=0 ){
    lines(xd, yd, col=norm.col)
  }
  
  return(invisible())
}##function plot.density.mcmcSTmodel
