#' Weighted NPMLE ofLongitudinal Gaussian Mean and Variances Model
#' 
#' A Kiefer-Wolfowitz procedure for ML estimation of a Gaussian model with
#' dependent mean and variance components and weighted longitudinal data.
#' This version assumes a general bivariate distribution for the mixing
#' distribution. The defaults use a rather coarse bivariate gridding.
#' 
#' @param y A vector of observations
#' @param id A strata indicator vector of the same length
#' @param w A vector of weights
#' @param u A vector of bin boundaries for the mean effects
#' @param v A vector of bin boundaries for the variance effects
#' @param ... optional parameters to be passed to KWDual to control optimization
#' @return A list consisting of the following components: 
#' 	\item{u}{midpoints of mean bin boundaries} 
#' 	\item{v}{midpoints of variance bin boundaries} 
#' 	\item{fuv}{the function values of the mixing density.} 
#' 	\item{logLik}{log likelihood value for mean problem} 
#' 	\item{du}{Bayes rule estimate of the mixing density means.} 
#' 	\item{dv}{Bayes rule estimate of the mixing density variances.} 
#' 	\item{status}{Mosek convergence status}
#' @author R. Koenker and J. Gu
#' @references Gu, J. and R. Koenker (2014) Heterogeneous Income Dynamics: An
#' Empirical Bayes Perspective, \emph{JBES}, forthcoming.
#' @seealso WTLVmix for an implementation assuming independent heterogeneity
#' @keywords nonparametric
#' @export
WGLVmix <- function(y, id, w, u = 30, v = 30, ...){

    n <- length(y)
    eps <- 1e-4
    if(missing(w)) w <- rep(1,n)
    wsum <- tapply(w, id, "sum")
    t <- tapply(w * y, id, "sum")/wsum
    m <- tapply(y, id, "length")
    r <- (m - 1)/2
    s <- (tapply(w * y^2, id, "sum") - t^2 * wsum)/(m - 1)
    n <- length(s)
    if(length(u) == 1) u <- seq(min(t) - eps, max(t) + eps, length = u)
    if(length(v) == 1) v <- seq(min(s) - eps, max(s) + eps, length = v)
    pu <- length(u)
    du <- diff(u)
    u <- (u[-1] + u[-pu])/2
    pu <- length(u)
    wu <- rep(1,n)/n
    pv <- length(v)
    dv <- diff(v)
    v <- (v[-1] + v[-pv])/2
    pv <- length(v)
    
    R <- outer(r*s,v,"/")  
    G <- outer(s * gamma(r),rep(1,pv))
    r <- outer((m - 1)/2, rep(1,pv))
    Av <- outer((exp(-R) * R^r)/G, rep(1,pu))
    Au <- dnorm(outer(outer(t, u, "-") * outer(sqrt(wsum),rep(1,pu)), sqrt(v), "/"))
    Au <- Au/outer(outer(1/sqrt(wsum),rep(1,pu)),sqrt(v))
    Au <- aperm(Au,c(1,3,2)) # permute Au indices to align with those of Av
    A <- Av * Au
    
    B <- NULL
    for (j in 1:pu) B <- cbind(B,A[,,j])
    duv = as.vector(kronecker(du, dv))
    uv <- expand.grid(theta = v, alpha = u)
    f <- KWDual(B, duv, wu, ...)
    fuv <- f$f
    g = f$g
    status <- f$status
    r <- (m-1)/2
    logK <- log(gamma(r)) - r * log(r) - 0.5 * log(wsum) - 
	r * log(2*pi) - log(s^(r-1)) + 0.5 * tapply(log(w), id, "sum")
    logLik <- sum(log(g)) + sum(logK)
    dx <- B%*%(uv[,2] * duv * fuv)/g  #Bayes rule for u: E(u|t,s)
    dy <- B %*% (uv[,1] * duv * fuv)/g  # Bayes rule for v: E(v|t,s)
    list(u = u, v = v, fuv = fuv, logLik = logLik, dx = dx, dy = dy,
	status = status)
}
