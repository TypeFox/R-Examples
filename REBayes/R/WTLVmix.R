#' NPMLE for Longitudinal Gaussian Means and Variances Model
#' 
#' A Kiefer-Wolfowitz NPMLE procedure for estimation of a Gaussian model with
#' independent mean and variance components with weighted longitudinal data.
#' This version exploits a Student t decomposition of the likelihood.
#' 
#' @param y A vector of observations
#' @param id A strata indicator vector indicating grouping of y
#' @param w A vector of weights corresponding to y
#' @param u A vector of bin boundaries for the mean effects
#' @param v A vector of bin boundaries for the variance effects
#' @param ... optional parameters to be passed to KWDual to control optimization
#' @return A list consisting of the following components: 
#' 	\item{u}{midpoints of the mean bin boundaries} 
#' 	\item{fu}{the function values of the mixing density of the means } 
#' 	\item{v}{midpoints of the variance bin boundaries} 
#' 	\item{fv}{the function values of the mixing density of the variances.} 
#' 	\item{logLik}{log likelihood value for mean problem} 
#' 	\item{du}{Bayes rule estimate of the mixing density means.} 
#' 	\item{dv}{Bayes rule estimate of the mixing density variances.} 
#' 	\item{status}{Mosek convergence status}
#' @author J. Gu and R. Koenker
#' @seealso WGLVmix for a more general bivariate mixing distribution version
#' @references Gu, J. and R. Koenker (2015) Empirical Bayesball Remixed, preprint
#' @keywords nonparametric
#' @export
WTLVmix <- function(y, id, w, u = 300, v = 300, ...) {

    n <- length(y)
    if(missing(w)) w <- rep(1,n)
    wsum <- tapply(w,id,"sum")
    t <- tapply(w*y,id,"sum")/wsum
    m <- tapply(y,id,"length")
    r <- (m-1)/2
    s <- (tapply(w*y^2,id,"sum") - t^2*wsum)/(m-1)
    n <- length(s)
    eps <- 1e-4
    if(length(u) == 1) 
	u <- seq(min(t) - eps, max(t) + eps, length = u)
    if(length(v) == 1) 
	v <- seq(min(s) - eps, max(s) + eps, length = v)
    pu <- length(u)    
    pv <- length(v)    
    du <- diff(u)
    u <- (u[-1] + u[-pu])/2
    wu <- rep(1,n)/n
    pu <- length(u)    

    Au <- dt(outer(t, u, "-") * outer(sqrt(wsum/s), rep(1, pu)), df = m - 1)
    Au <- Au/outer(sqrt(s/wsum), rep(1, pu))
    f <- GVmix(s, m, v = v, ...)
    dv <- diff(v)
    v <- (v[-1] + v[-pv])/2
    pv <- length(v)    
    fv <- f$y
    status <- f$status
    f <- KWDual(Au, du, wu, ...)
    fu <- f$f
    status <- c(f$status, status)

    # Likelihood and Bayes rule computation
    R <- outer(r*s,v,"/")
    sgamma <- outer(s * gamma(r),rep(1,pv))
    r <- outer((m - 1)/2, rep(1,pv))
    Av <- outer((exp(-R) * R^r)/sgamma, rep(1,pu))
    Au <- dnorm(outer(outer(t, u, "-") * outer(sqrt(wsum),rep(1,pu)), sqrt(v), "/"))
    Au <- Au/outer(outer(1/sqrt(wsum),rep(1,pu)),sqrt(v))
    Au <- aperm(Au,c(1,3,2)) # permute Au indices to align with those of Av
    A <- Av * Au #dim: n*pv*pu
    au<-matrix(0,n,pu)  # marginalize out v
    for (i in 1:pu) 
	au[,i]<-A[,,i]%*%(dv*fv)
    av <- matrix(0, n, pv)  #marginalize out u
    for (i in 1:pv)
        av[,i] <- A[,i,] %*%(du * fu)
    g <- au%*%(du*fu)  # marginal density for vector of y_i
    dx <- au%*%(u * du * fu)/g  # Bayes rule E(u|t,s)
    dy <- av %*%(v * dv * fv)/g # Bayes rule E(v|t,s)
    # an extra factor in the likelihood
    r <- (m-1)/2
    logK <- log(gamma(r)) - r*log(r) - 0.5 * log(wsum) - 
	r*log(2*pi) - log(s^(r-1)) + 0.5 * tapply(log(w), id, "sum")
    logLik <- sum(logK) + sum(log(g))
    list(u = u, fu = fu, v = v, fv = fv, logLik = logLik, dx = dx, dy = dy,
	status = status)
}

