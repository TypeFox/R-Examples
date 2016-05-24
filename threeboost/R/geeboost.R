#' @title GEEBoost
#' 
#' @description
#' Thresholded boosting for correlated data via GEE
#' 
#' @export
#' 
#' @details
#' This function implements thresholded EEBoost for the Generalized Estimating Equations.
#' The arguments are consistent with those used by \code{geepack}.
#' 
#' @param Y Vector of (presumably correlated) outcomes
#' @param X Matrix of predictors
#' @param id Index indicating clusters of correlated observations
#' @param family Outcome distribution to be used. \code{"gaussian"} (the default), \code{"binomial"}, and \code{"poisson"} are currently implemented.
#' @param corstr Working correlation structure to use. \code{"ind"} (for independence, the default) and \code{"exch"} (for exchanageable) are currently implemented.
#' @param traceplot Option of whether or not to produce a traceplot of the coefficient values. See \code{\link{coef_traceplot}} for details.
#' @param ... Additional arguments to be passed to the \code{\link{threeboost}} function. See the \code{\link{threeboost}} help page for details. 
#' 
#' @return A list with three entries:
#' \itemize{
#'  \item \code{coefmat} A matrix with \code{maxit} rows and \code{ncol(X)} columns, 
#'  with each row containing the parameter vector from an iteration of EEBoost.
#'  \item \code{QICs} A vector of QICs computed from the coefficients.
#'  \item \code{final.model} The coefficients corresponding to the model (set of coefficients)
#'  yielding the smallest QIC.
#'  }
#' 
#' @examples
#' # Generate some test data
#' library(mvtnorm)
#' library(Matrix)
#' n <- 30
#' n.var <- 50
#' clust.size <- 4
#' 
#' B <- c(rep(2,5),rep(0.2,5),rep(0.05,10),rep(0,n.var-20))
#' mn.X <- rep(0,n.var)
#' sd.X <- 0.5
#' rho.X <- 0.3
#' cov.sig.X <- sd.X^2*((1-rho.X)*diag(rep(1,10)) + rho.X*matrix(data=1,nrow=10,ncol=10))
#' sig.X <- as.matrix( Matrix::bdiag(lapply(1:(n.var/10),function(x) { cov.sig.X } ) ) )
#' sd.Y <- 0.5
#' rho.Y <- 0.3
#' indiv.Sig <- sd.Y^2*( (1-rho.Y)*diag(rep(1,4)) + rho.Y*matrix(data=1,nrow=4,ncol=4) )
#' sig.list <- list(length=n)
#' for(i in 1:n) { sig.list[[i]] <- indiv.Sig }
#' Sig <- Matrix::bdiag(sig.list)
#' indiv.index <- rep(1:n,each=clust.size)
#' sig.Y <- as.matrix(Sig)
#' 
#' if(require(mvtnorm)) {
#' X <- mvtnorm::rmvnorm(n*clust.size,mean=mn.X,sigma=sig.X)
#' mn.Y <- X %*% B
#' Y <- mvtnorm::rmvnorm(1,mean=mn.Y,sigma=sig.Y) ## Correlated continuous outcomes
#' expit <- function(x) { exp(x) / (1 + exp(x)) }
#' ## Correlated binary outcomes
#' Y.bin <- rbinom(n*clust.size,1,p=expit(mvtnorm::rmvnorm(1,mean=mn.Y,sigma=sig.Y))) 
#' Y.pois <- rpois(length(Y),lambda=exp(mn.Y)) ## Correlated Poisson outcomes
#' } else { stop('Need mvtnorm package to generate correlated data.')}
#' 
#' ## Run EEBoost (w/ indep working correlation)
#' results.lin <- geeboost(Y,X,id=indiv.index,maxit=1000)
#' \dontrun{
#'  results.bin <- geeboost(Y.bin,X,id=indiv.index,family="binomial",maxit=1000)
#'  results.pois <- geeboost(Y.pois,X,id=indiv.index,family="poisson",maxit=1000,traceplot=TRUE)
#' }
#' 
#' print(results.lin$final.model)
#' 
#' @seealso \code{\link{threeboost}}
#' 
#' Wolfson, J. \href{http://pubs.amstat.org/doi/abs/10.1198/jasa.2011.tm10098}{EEBoost: A general method for prediction and variable selection using estimating equations.} Journal of the American Statistical Association, 2011.
geeboost <- function(Y,X,id=1:length(Y),family="gaussian",corstr="ind",traceplot=FALSE,...) {
  
  main.fn <- switch(family,
                    gaussian = ee.GEELin,
                    binomial = ee.GEEBin,
                    poisson = ee.GEEPois)
  aux.fn <- switch(family,
                   gaussian = ee.GEELin.aux,
                   binomial = ee.GEEBin.aux,
                   poisson = ee.GEEPois.aux)
  mu.fn <- switch(family,
                  gaussian = mu.Lin,
                  binomial = mu.Bin,
                  poisson = mu.Pois)
  g.fn <- function(family,
                   gaussian = g.Lin,
                   binomial = g.Bin,
                   poisson = g.Pois)
  v.fn <- switch(family,
                 gaussian = v.Lin,
                 binomial = v.Bin,
                 poisson = v.Pois)
  EE.fn <- switch(corstr,
                  ind = function(Y,X,b) {
                    main.fn(as.vector(Y),X,b,aux=aux.fn,id=id,corstr="ind") },
                  exch = function(Y,X,b) {
                    main.fn(as.vector(Y),X,b,aux=aux.fn,id=id,corstr="exch") } )
  
  QIC.fn <- switch(family,
                   gaussian=QIC.lin,
                   binomial=QIC.bin,
                   poisson=QIC.pois)
  ##boost.out <- threeboost(Y=as.vector(Y),X=X,EE.fn=EE.fn)
  boost.out <- threeboost(Y=as.vector(Y),X=X,EE.fn=EE.fn,...)
  intcpts <- apply(boost.out,1,function(b) {
    aux.par <- aux.fn(Y=as.vector(Y),X=X,b=b,id=id)
    return(aux.par[3])
  })
  coefmat <- cbind(intcpts,boost.out)

  X.int <- cbind(rep(1,nrow(X)),X)
  QICs <- apply(coefmat,1,function(b) {
    QIC(as.vector(Y),X.int,b,family=family)
  })
  final.model <- coefmat[which.min(QICs),]
  
  if(traceplot==TRUE) { coef_traceplot(coefmat[,-1]) } ## Exclude the intercept in the traceplot
  
  return(list(coefmat=coefmat,QICs=QICs,final.model=final.model))
} 
