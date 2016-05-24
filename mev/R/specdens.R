
#' Rank-based transformation to angular measure
#'
#' The method uses the pseudo-polar transformation for suitable norms, transforming
#' the data to pseudo-observations, than marginally to unit Frechet or unit Pareto.
#' Empirical or Euclidean weights are computed and return alongside of the angular and
#' radial sample
#'
#' @param x an \code{n} by \code{d} sample matrix
#' @param Rnorm the norm for the radial component
#' @param Anorm the norm for the angular component. \code{arctan} is only implemented for \eqn{d=2}
#' @param marg choice of marginal transformation, either to Frechet or Pareto scale
#' @param wgt weighting function for the equation. Can be based on Euclidean or empirical likelihood for the mean
#' @return a list with arguments \code{ang} for the \eqn{d-1} pseudo-angular sample, \code{rad} with the radial component
#' and \code{wts} if \code{Rnorm} is set to "\code{l1}" (default).
#' @author Leo Belzile
#' @references Einmahl, J.H.J. and J. Segers (2009). Maximum empirical likelihood estimation of the spectral measure of an extreme-value distribution, \emph{Annals of Statistics}, \bold{37}(5B), 2953--2989.
#' @references de Carvalho, M. and B. Oumow and J. Segers and M. Warchol (2013). A Euclidean likelihood estimator for bivariate tail dependence, \emph{Comm. Statist. Theory Methods}, \bold{42}(7), 1176--1192.
#' @references Owen, A.B. (2001). \emph{Empirical Likelihood}, CRC Press, 304p.
#' @export
#' @return a list with components
#' \itemize{
#' \item \code{ang} matrix of pseudo-angular observations
#' \item \code{rad} vector of radial contributions
#' \item \code{wts} empirical or Euclidean likelihood weights for observations
#' }
#'
#' @examples
#' x <- rmev(n=25, d=3, param=0.5, model="log")
#' wts <- angmeas(x=x, Rnorm="l1", Anorm="l1", marg="Frechet", wgt="Empirical")
#' wts2 <- angmeas(x=x, Rnorm="l2", Anorm="l2", marg="Pareto", wgt="Euclidean")
angmeas <- function(x, Rnorm=c("l1","l2","linf"), Anorm=c("l1","l2","linf","arctan"),
  marg=c("Frechet","Pareto"), wgt=c("Euclidean","Empirical")){
  if (!is.matrix(x)) x <- rbind(x, deparse.level = 0L)
    marg <- match.arg(marg)
    S <- switch(marg,
      Frechet=-1/log(apply(x, 2, rank, na.last = "keep", ties.method = "random")/(nrow(x) +  1)),
      Pareto = 1/(1-apply(x, 2, rank, na.last = "keep", ties.method = "random")/(nrow(x) +  1))
    )
    R <- switch(Rnorm,
      l1   = rowSums(S),
      l2   = apply(S,1,function(x){sqrt(sum(x^2))}),
      linf = apply(S,1,max)
    )
    #Ordering is same regardless of above norm
    ordR <- order(R)
    S <- S[ordR,]
    R <- R[ordR]
    if(Rnorm==Anorm){
    ang <- S[,-ncol(S)]/R
    } else if(Anorm=="arctan"){
      if(ncol(S)!=2){
        stop("Invalid norm for sample, arctan transformation only for bivariate samples")
        }
     ang <- atan(S[,2]/S[,1])
    } else{
     ang <- switch(Anorm,
        l1   = S[,-ncol(S)]/rowSums(S),
        l2   = S[,-ncol(S)]/apply(S,1,function(x){sqrt(sum(x^2))}),
        linf = S[,-ncol(S)]/apply(S,1,max)
     )
    }
    rm(ordR, S)
      if(wgt=="Euclidean"){
        return(list(ang=ang,rad=R, wts=as.vector(.EuclideanWeights(ang,rep(1/(ncol(ang)+1), ncol(ang))))))
      } else if(wgt=="Empirical"){
      		scel.fit <- .emplik(z=ang,mu=rep(1/(ncol(ang)+1), ncol(ang)), lam=rep(0,ncol(ang)), eps=1/nrow(ang))
      		if(scel.fit$conv){
        return(list(ang=ang,rad=R, wts=as.vector(scel.fit$wts)))
      		} else{
      			warning("Self-concordant empirical likelihood for the mean did not converge.
											Verify arguments or else try increasing the number of iterations.")
      	return(list(ang=ang,rad=R))
      		}
      }
}

#' Self-concordant empirical likelihood for a vector mean
#'
#' @param dat \code{n} by \code{d} matrix of \code{d}-variate observations
#' @param mu  \code{d} vector of hypothesized mean of \code{dat}
#' @param lam  starting values for Lagrange multiplier vector, default to zero vector
#' @param eps  lower cutoff for \eqn{-\log}{-log}, with default \code{1/nrow(dat)}
#' @param M upper cutoff for \eqn{-\log}{-log}.
#' @param thresh  convergence threshold for log likelihood (default of \code{1e-30} is agressive)
#' @param itermax  upper bound on number of Newton steps.
#' @export
#' @author Art Owen, \code{C++} port by Leo Belzile
#' @references Owen, A.B. (2013). Self-concordance for empirical likelihood, \emph{Canadian Journal of Statistics}, \bold{41}(3), 387--397.
#' @return a list with components
#' #' \itemize{
#'  \item \code{logelr} log empirical likelihood ratio.
#'  \item \code{lam} Lagrange multiplier (vector of length \code{d}).
#'  \item \code{wts} \code{n} vector of observation weights (probabilities).
#'  \item \code{conv} boolean indicating convergence.
#'  \item \code{niter} number of iteration until convergence.
#'  \item \code{ndec} Newton decrement.
#'  \item \code{gradnorm} norm of gradient of log empirical likelihood.
#' }
emplik <- function(dat, mu=rep(0, ncol(dat)), lam = rep(0, ncol(dat)), eps = 1/nrow(dat), M=1e30, thresh=1e-30, itermax=100){
	if(is.infinite(M))		M = 1e30;
	.emplik(z=dat, mu=mu, lam, eps=eps, M=M, thresh=thresh, itermax=itermax)
}
