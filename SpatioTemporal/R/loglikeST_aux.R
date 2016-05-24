#######################################################
## HELPER FUNCTIONS FOR LIKELIHOOD, SIZE, NAMES, ETC ##
#######################################################
##Functions in this file:
## loglikeSTdim     EX:ok
## loglikeSTnames   EX:ok
## loglikeSTgetPars EX:ok

##' Function that computes the dimension of several objects in a \code{STmodel}
##' object.
##'
##' @title Dimensions of the STmodel Structure
##' @param STmodel \code{STmodel} object for which dimensions are to
##'   be computed.
##'
##' @return list containing:
##'   \item{T}{Number of observation times.}
##'   \item{m}{Number of temporal basis functions, including the intercept.}
##'   \item{n}{Number of distinct locations in the data.}
##'   \item{n.obs}{Number of observed locations.}
##'   \item{p}{vector of length \code{m}; number of geographic covariates
##'            for each temporal basis functions.}
##'   \item{L}{Number of spatio-temporal covariates}
##'   \item{npars.beta.covf}{vector of length \code{m}; number of parameters for each
##'                          covariance-function for the beta-fields.}
##'   \item{npars.beta.tot}{vector of length \code{m}; total number of
##'                         parameters for each 
##'                         beta-fields, including nugget(s).}
##'   \item{npars.nu.covf,npars.nu.tot}{number of parameters for the nu-field, same
##'                                     distinction as above.}
##'   \item{nparam}{Total number of parameters, including regression parameters.}
##'   \item{nparam.cov}{Number of covariance parameters.}
##' 
##' @examples 
##' ##load the data
##' data(mesa.model)
##' 
##' ##compute dimensions for the data structure
##' loglikeSTdim(mesa.model)
##'
##' @author Johan Lindström
##' 
##' @family STmodel functions
##' @family likelihood utility functions
##' @export
loglikeSTdim <- function(STmodel){
  ##check class belonging
  stCheckClass(STmodel, "STmodel", name="STmodel")
  
  out <- list()
  ##nbr of time points
  out$T <- length( STmodel$nt )
  ##nbr of temporal basis fnct. incld. the constant
  out$m <- dim(STmodel$F)[2]
  ##nbr of sites
  out$n <- dim(STmodel$LUR.all[[1]])[1]
  ##nbr of observed sites
  out$n.obs <- dim(STmodel$LUR[[1]])[1]

  ##nbr land use regressions for each temporal basis fnct.
  out$p <-  sapply(STmodel$LUR, function(x){dim(x)[2]})

  ##number of spatio-temporal covariates
  out$L <- length(STmodel$ST.list)

  ##number of parameters for the covariance functions
  out$npars.beta.covf <- sapply( parsCovFuns(STmodel$cov.beta$covf, list=TRUE),
                                length)
  out$npars.beta.tot <- (out$npars.beta.covf +  STmodel$cov.beta$nugget)
  out$npars.nu.covf <- length( parsCovFuns(STmodel$cov.nu$covf, list=FALSE) )
  out$npars.nu.tot <- out$npars.nu.covf + STmodel$cov.nu$random.effect
  out$npars.nu.tot <- out$npars.nu.tot + dim(STmodel$cov.nu$nugget.matrix)[2]
  ##number of unknown parameters
  out$nparam.cov <- sum(out$npars.beta.tot)+out$npars.nu.tot
  out$nparam <- out$nparam.cov + sum(out$p) + out$L
  ##return
  return( out )
}##function loglikeSTdim


##' Function that creates a character vector with names for the parameters
##' expected by log-likelihood functions. Names are created by extracting names
##' from the \code{STmodel} structure.
##'
##' @title Create Names for Log-likelihood Parameters for STmodel objects
##' @param STmodel \code{STmodel} object for which parmeter names are to
##'   be computed.
##' @param all compute all parameter names (regression and covariance) or
##'   only covariance parameters.
##' @return Returns names of the parameters expected by the log-likelihood
##'   functions. Regerssion parameter names start with gamma/alpha
##'   (spatio-temporal/geographic covariate), followed by name of beta-field,
##'   and the name of covariate. The covariance parameters follow,
##'   log (reminder that parameter is log-scale), covariance parameter name,
##'   name of field, type of covariance function.
##' 
##' @examples 
##' ##load the data
##' data(mesa.model)
##' 
##' ##Find out in which order parameters should be given
##' loglikeSTnames(mesa.model)
##' ##...and for only the covariance parameters.
##' loglikeSTnames(mesa.model, FALSE)
##'
##' @author Johan Lindström
##' 
##' @family likelihood utility functions
##' @export
loglikeSTnames <- function(STmodel, all=TRUE){
  ##check class belonging
  stCheckClass(STmodel, "STmodel", name="STmodel")
  
  dimensions <- loglikeSTdim(STmodel)
  ##create a blank vector of the right length
  out <- NULL
  ##add appropriate names to said vector
  if(all){
    if(dimensions$L!=0){ ##regression coefficients for the model
      out <- c(out, paste("gamma", colnames(STmodel$ST), sep="."))
    }
    tmp <- lapply(STmodel$LUR, colnames)
    for(i in 1:length(tmp)){
      tmp[[i]] <- paste("alpha",names(STmodel$LUR)[i],tmp[[i]],sep=".")
    }
    out <- c(out, unlist(tmp))
  }##if(all)

  ##covariance for beta
  tmp <- parsCovFuns(STmodel$cov.beta$covf, list=TRUE)
  for(i in 1:length(tmp)){
    tmp[[i]] <- c(tmp[[i]], switch(STmodel$cov.beta$nugget[i], "nugget"))
    tmp[[i]] <- paste("log", tmp[[i]], names(STmodel$LUR)[i], names(tmp)[i],
                      sep=".")
  }
  out <- c(out, unlist(tmp))

  ##covariance for nu
  if( dim(STmodel$cov.nu$nugget.matrix)[2]==0 ){
    tmp <- c(parsCovFuns(STmodel$cov.nu$covf, list=FALSE),
             switch(STmodel$cov.nu$random.effect, "random.effect"))
  }else{
    tmp <- c(parsCovFuns(STmodel$cov.nu$covf, list=FALSE),
             paste("nugget", colnames(STmodel$cov.nu$nugget.matrix), sep="."),
             switch(STmodel$cov.nu$random.effect, "random.effect"))
  }
  tmp <- paste("nu", "log", tmp, STmodel$cov.nu$covf, sep=".")

  ##concatonate, remove names and return
  out <- c(out, tmp)
  names(out) <- NULL
  return(out)
}##function loglikeSTnames

##' Extracts parameters for the log-likelihood from a parameter vector and
##' separates regression parameters and \emph{log}-covariance parameters.
##'
##' @title Extract Parameters from a Vector
##' @param x A vector containing regression (optionally)  and
##'   \emph{log}-covariance parameters. The ordering of
##'   \strong{has to be exactly} that indicated by
##'   \code{\link{loglikeSTnames}}.
##' @param STmodel STmodel \code{STmodel} object describing the problem.
##' 
##' @return list containing:
##'   \item{gamma}{Regression coefficients for the spatio-temporal covariate(s).}
##'   \item{alpha}{A list of regression coefficients for geographic covariates.}
##'   \item{cov.beta}{A list containg a lists of pars and vector of nuggets.
##'                   See \code{\link{makeSigmaB}}.}
##'   \item{cov.nu}{A list of covariance parameters for the nu-field, as
##'     \code{pars}, \code{nugget} and \code{random.effect} respectively.}
##' Covariance parameters are also back-transformed from log-scale.
##'
##' @example Rd_examples/Ex_loglikeSTgetPars.R
##' 
##' @author Johan Lindström
##' 
##' @family likelihood utility functions
##' @export
loglikeSTgetPars <- function(x, STmodel){
  ##check class belonging
  stCheckClass(STmodel, "STmodel", name="STmodel")
  ##compute dimensions
  dimensions <- loglikeSTdim(STmodel)

  if( length(x)!=dimensions$nparam && length(x)!=dimensions$nparam.cov ){
    stop( paste("Parameter vector 'x' is", length(x), "should be",
                dimensions$nparam.cov, "or", dimensions$nparam) )
  }
  Ind <- 0
  if( length(x)==dimensions$nparam ){
    ##regression coefficients for model output
    if( dimensions$L!=0 ){ 
      gamma <- matrix( x[Ind + (1:dimensions$L)], ncol=1)
      Ind <- Ind + dimensions$L
    }else{
      gamma <- double(0)
    }
    ##land use regression coefficients for each of the temporal basis functions
    alpha <- vector("list", dimensions$m)
    for( i in 1:dimensions$m){
      alpha[[i]] <- matrix( x[Ind + (1:dimensions$p[i])], ncol=1)
      Ind <- Ind + dimensions$p[i]
    }
  }else{
    alpha <- double(0)
    gamma <- double(0)
  }
  ##covariance for beta (and invert log trans form, i.e. exp)
  cov.beta <- list(pars=vector("list", dimensions$m),
                   nugget=double(dimensions$m))
  for(i in 1:dimensions$m){
    ##covariance parameters
    if( dimensions$npars.beta.covf[i]!=0 ){
      cov.beta$pars[[i]] <- exp( x[Ind+(1:dimensions$npars.beta.covf[i])] )
      Ind <- Ind+dimensions$npars.beta.covf[i]
    }else{
      cov.beta$pars[[i]] <- numeric(0)
    }
    ##nugget
    if( STmodel$cov.beta$nugget[i] ){
      cov.beta$nugget[i] <- exp( x[Ind+1] )
      Ind <- Ind+1
    }else{
      cov.beta$nugget[i] <- 0
    }
  }##for(i in 1:dimensions$m)

  ##covariance for nu
  cov.nu <- list()
  ##covariance parameters
  if( dimensions$npars.nu.covf!=0 ){
    cov.nu$pars <- exp( x[Ind+(1:dimensions$npars.nu.covf)] )
    Ind <- Ind+dimensions$npars.nu.covf
  }else{
    cov.nu$pars <- numeric(0)
  }
  ##nugget
  if( dim(STmodel$cov.nu$nugget.matrix)[2]==0 ){
    cov.nu$nugget <- matrix(0, dim(STmodel$cov.nu$nugget.matrix)[1], 1)
  }else{
    tmp <- matrix(x[Ind+(1:dim(STmodel$cov.nu$nugget.matrix)[2])], ncol=1)
    cov.nu$nugget <- exp( STmodel$cov.nu$nugget.matrix %*% tmp )
    Ind <- Ind+dim(STmodel$cov.nu$nugget.matrix)[2]
  }
  rownames(cov.nu$nugget) <- rownames(STmodel$cov.nu$nugget.matrix)
  ##random effect
  if( STmodel$cov.nu$random.effect ){
    cov.nu$random.effect <- exp( x[Ind+1] )
  }else{
    cov.nu$random.effect <- 0
  }

  return( list(gamma=gamma, alpha=alpha, cov.beta=cov.beta, cov.nu=cov.nu) )
}##function loglikeSTgetPars

