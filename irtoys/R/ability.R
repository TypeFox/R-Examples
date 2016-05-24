# 3PL log posterior for one person / response vector
# r=responses, p=parm list, x=ability, mu, sigma
llf = function(x,r,p,mu,sigma,method) {
	pr = p[,3] + (1.0 - p[,3])/(1.0 + exp(p[,1]*(p[,2] - x)))
	pr = pmax(pr, .00001); pr = pmin(pr, .99999)
	ll = r*log(pr) + (1-r)*log(1.0-pr)
	lf = sum(ll)
	if (method != "ML") lf = lf + log(dnorm(x,mu,sigma)) 
  return(lf)
} 

mle.one = function(resp, ip, mu=mu, sigma=sigma, method=method) {                                                            
    cc = !is.na(resp)                                        
    resp = resp[cc]                                          
    ip = ip[cc, , drop=FALSE]                                             
    n = length(resp)                                         
    if (n < 1) return(c(NA, NA, 0))                                 
    est = optimize(llf, lower = -4, upper = 4, maximum = TRUE, 
        r = resp, p = ip, mu = mu, sigma = sigma, method = method)$maximum
    ti = tif(ip, est)$f
    if (method != "ML") ti = ti + 1/(sigma * sigma)
    sem = sqrt(1/ti)
    return(c(est, sem, n))
}


#' Normal quadrature points and weights
#' 
#' Quadrature points and weights based on the Normal distribution. Quadrature
#' objects are used when estimating abilities with \code{eap} and for some of
#' the scaling methods in \code{sca}.
#' 
#' 
#' @param n Number of quadrature points
#' @param lower Lower boundary
#' @param upper Upper boundary
#' @param mu Mean
#' @param sigma Standard deviation
#' @param scaling Determines the way in which non-default values of \code{mu}
#' and \code{sigma} are applied. When \code{scaling="points"}, the quadrature
#' points are rescaled, otherwise the quadrature weights are adapted.
#' @return A list of: \item{quad.points}{A vector of \code{n} quadrature
#' points} \item{quad.weights}{A vector of the corresponding quadrature
#' weights}
#' @author Ivailo Partchev
#' @export
#' @seealso \code{\link{read.qu.icl}}, \code{\link{eap}}, \code{\link{sca}}
#' @keywords models
#' @examples
#' 
#' quad <- normal.qu(n=20)
#' 
normal.qu = function(n=15,lower=-4,upper=4,mu=0,sigma=1,scaling="points"){
  if (upper<=lower || sigma<=0 || n<3) stop("bad argument")
  qp=seq(lower,upper,length.out=n)
  if(scaling=="points") {
  	qw=dnorm(qp,0,1)
  	qw=qw/sum(qw)
  	qp=qp*sigma+mu
  } else {
  	qw=dnorm(qp,mu,sigma)
  	qw=qw/sum(qw)
  }
  return(list(quad.points=qp, quad.weights=qw))
}


#' Maximum likelihood and Bayes Modal estimation of ability
#' 
#' Estimates the value of the latent variable ("ability") for each person by
#' direct optimization
#' 
#' 
#' @param resp A matrix of responses: persons as rows, items as columns,
#' entries are either 0 or 1, no missing data
#' @param ip Item parameters: a matrix with one row per item, and three
#' columns: [,1] item discrimination \eqn{a}, [,2] item difficulty \eqn{b}, and
#' [,3] asymptote \eqn{c}.
#' @param mu Mean of the apriori distribution. Ignored when \code{method="ML"}.
#' Default is 0.
#' @param sigma Standard deviation of the apriori distribution. Ignored when
#' \code{method="ML"}. Default is 1.
#' @param method \code{"ML"} for maximum likelihood or \code{"BM"} for Bayes
#' Modal estimation. Default is \code{"ML"}, in which case any values for
#' \code{mu} and \code{sigma} will be ignored.
#' @return A matrix with the ability estimates in column 1 and their standard
#' errors of measurement (SEM) in column 2, and the number of non-missing
#' responses in column 3
#' @author Ivailo Partchev
#' @seealso \code{\link{eap}}
#' @keywords models
#' @export
#' @examples
#' 
#' th.mle <- mlebme(resp=Scored, ip=Scored2pl$est)
#' 
mlebme = function(resp, ip, mu=0, sigma=1, method="ML") {
 if (is.null(dim(resp))) dim(resp) = c(1,length(resp))
 if (is.null(dim(ip))) stop("item parameters not a matrix")
 if (nrow(ip) != ncol(resp)) stop("responses - item parameters mismatch")
 np = nrow(resp)
 o = sapply(1:np, function(i) mle.one(resp=resp[i,], 
    ip=ip, mu=mu, sigma=sigma, method=method))
 rownames(o) = c("est","sem","n")
 return(t(o)) 
}
# bias-corrected (Warm's) estimate for one person
bce.one = function(resp, ip) {                                                            
    cc = !is.na(resp)                                        
    resp = resp[cc]                                          
    ip = ip[cc, , drop=FALSE]                                             
    n = length(resp)                                         
    if (n < 1) return(c(NA, NA, 0))                                 
    est = uniroot(scf, re=resp, p=ip, lower=-10, upper=10)$root		
    ev = bcv(est, resp, ip)
	return(c(est, sqrt(ev), n))
}

# variance of the Warm estimator
bcv = function(x,r,p) {
  i = iif(p, x)$f
  p[,3] = 0
  q = irf(p, x)$f
  isum = sum(i)
  jsum = sum(i * p[, 1] * (1 - 2 * q))
  return(1/isum + jsum^2/(4 * isum^4))
}
	
# score function (Warm's estimates for the 3PL 
# are unwieldy for direct optimization so use 1st deriv)
# r=responses, p=parm list, x=ability
scf = function(x,re,p) {
	three = any(p[,3] > 0)
	lgt = exp(p[,1] * (x - p[,2]))
	pr = lgt / (1 + lgt)
	z = re - pr
	if (three) z = z - p[,3]*re / (p[,3] + lgt)
	sm = sum(p[,1]*z)
	if (three) {
		pr3 = p[,3] + (1 - p[,3])*pr
		ii = p[,1]^2 / pr3 * (1 - pr3) * pr^2
	} else {
		ii = p[,1]^2 * pr * (1 - pr)
	}
	isum = sum(ii)
	jsum = sum(ii * p[,1] * (1 - 2*pr))
	return(sm + jsum / (isum*2))
}

# Bias-corrected (aka Warm's) ability estimates


#' Bias-corrected (Warm's) estimates of ability
#' 
#' Weighted likelihood estimates (WLE) of ability, designed to remove the first
#' order bias term from the ML estimates. WLE are finite for response patterns
#' consisting of either uniformly wrong or uniformly correct responses.
#' 
#' 
#' @param resp A matrix of responses: persons as rows, items as columns,
#' entries are either 0 or 1, no missing data
#' @param ip Item parameters: a matrix with one row per item, and three
#' columns: [,1] item discrimination \eqn{a}, [,2] item difficulty \eqn{b}, and
#' [,3] asymptote \eqn{c}.
#' @return A matrix with the ability estimates in column 1, and their standard
#' errors of measurement (SEM) in column 2, and the number of non-missing
#' reponses in column 3
#' @author Ivailo Partchev
#' @seealso \code{\link{mlebme}}, \code{\link{eap}}
#' @references Warm T.A. (1989) Weighted Likelihood Estimation of Ability in
#' Item Response Theory. Psychometrika, 54, 427-450.
#' @keywords models
#' @export
#' @examples
#' 
#' th.bce <- wle(resp=Scored, ip=Scored2pl$est)
#' 
wle = function(resp, ip) {
 if (is.null(dim(resp))) dim(resp) = c(1,length(resp))
 if (is.null(dim(ip))) stop("item parameters not a matrix")
 if (nrow(ip) != ncol(resp)) stop("responses - item parameters mismatch")
 np = nrow(resp)
 o = sapply(1:np, function(i) bce.one(resp=resp[i,], ip=ip))
 rownames(o) = c("est","sem","n")
 return(t(o)) 
}
# 3PL EAP ability estimate for one person:
# r=responses, p=parm list, u=quad list)
eap.one = function(r, p, qp, qw) {
  cc = !is.na(r)
  r  = r[cc]
  p  = p[cc,,drop=FALSE]
  n  = length(r)
  if (n < 1) return(c(NA, NA, 0))
  ll = sapply(qp, llf, r=r, p=p, mu=NULL, sigma=NULL, method="ML")
  wl = exp(ll)*qw
  swl = sum(wl)
  x  = sum(wl*qp)/swl
  dev = qp - x
  sem = sqrt(sum(wl*dev*dev)/swl)
  return(c(x,sem,n))
}


#' EAP estimation of ability
#' 
#' Estimates the expectation of the posterior distribution of the latent
#' variable ("ability") for each person.
#' 
#' 
#' @param resp A matrix of responses: persons as rows, items as columns,
#' entries are either 0 or 1, no missing data
#' @param ip Item parameters: a matrix with one row per item, and three
#' columns: [,1] item discrimination \eqn{a}, [,2] item difficulty \eqn{b}, and
#' [,3] asymptote \eqn{c}.
#' @param qu A quadrature object produced with \code{\link{normal.qu}} or read
#' in with \code{\link{read.qu.icl}}
#' @return A matrix with the ability estimates in column 1, and their standard
#' errors of measurement (SEM) in column 2, and the number of non-missing
#' reponses in column 3
#' @author Ivailo Partchev
#' @seealso \code{\link{mlebme}}, \code{\link{normal.qu}},
#' \code{\link{read.qu.icl}}
#' @keywords models
#' @export
#' @examples
#' 
#' th.eap <- eap(resp=Scored, ip=Scored2pl$est, qu=normal.qu())
#' 
eap = function(resp, ip, qu) {
  if (is.null(dim(resp))) dim(resp) = c(1,length(resp))
  if (is.null(dim(ip))) stop("item parameters not a matrix")
  if (nrow(ip) != ncol(resp)) stop("responses - item parameters mismatch")
  np = nrow(resp)
  qp = qu$quad.points
  qw = qu$quad.weights
  o  = sapply(1:np, function(i) eap.one(r=resp[i,], p=ip, qp, qw),USE.NAMES=FALSE)
  rownames(o) = c("est","sem","n")
  return(t(o))
}

like = function(x, r, p, mu=0, s=1, log=FALSE, post=TRUE) {
  pr = irf(p,x)$f
	pr = pmax(pr, .00001); pr = pmin(pr, .99999)
  ll = log(pr) %*% r + log(1 - pr) %*% (1-r)
  if (post) 
    if (log) ll=ll+dnorm(x,mu,s,log=TRUE) else ll=exp(ll)*dnorm(x,mu,s)
  else if (!log) ll=exp(ll)
  return(ll)
}

ddf = function(x,r,p,d,mu,s) 
  log(like(x,r,p,mu=mu,s=s,post=TRUE)/d) - dt(x,df=3,log=TRUE) 

# rejection sampling for one person
dpv.one = function(resp, ip, n=5, mu, s) {
  cc = !is.na(resp)
  resp = resp[cc]
  ip   = ip[cc,]
  if (length(resp) < 1) return(rep(NA,n))
  d   = integrate(like, lower=-6, upper=6, p=ip, r=resp, mu=mu, s=s, post=TRUE)$value 
  dd  = optimize(f=ddf, c(-6,6), r=resp, p=ip, d=d, mu=mu, s=s, maximum=TRUE)$objective
  pv = rep(0,n)
  k  = 0
  repeat {
    th = rt(1, df=3)
    lf = log(like(th, r=resp, p=ip, mu=mu, s=s, post=TRUE) / d)  
    lg = dt(th, df=3, log=TRUE)
    prob = exp(lf - lg - dd)
    if (runif(1) < prob) {k = k+1; pv[k] = th}
    if (k==n) break
  }
  return(pv)
}

#' Draw plausible values
#' 
#' Draws (by rejection sampling) plausible values from the posterior
#' distribution of ability
#' 
#' 
#' @param resp A matrix of responses: persons as rows, items as columns,
#' entries are either 0 or 1, no missing data
#' @param ip Item parameters: a matrix with one row per item, and three
#' columns: [,1] item discrimination \eqn{a}, [,2] item difficulty \eqn{b}, and
#' [,3] asymptote \eqn{c}.
#' @param mu Mean of the apriori distribution. Ignored when \code{method="ML"}.
#' Default is 0.
#' @param sigma Standard deviation of the apriori distribution. Ignored when
#' \code{method="ML"}. Default is 1.
#' @param n The number of plausible values to draw for each person (default is
#' 5).
#' @return A matrix with \code{n} columns
#' @author Ivailo Partchev
#' @seealso \code{\link{mlebme}}, \code{\link{eap}}
#' @keywords models
#' @export
#' @examples
#' 
#' plval <- dpv(resp=Scored, ip=Scored2pl$est)
#' 
dpv = function(resp, ip, mu=0, sigma=1, n=5) {
 if (is.null(dim(resp))) dim(resp) = c(1,length(resp))
 if (is.null(dim(ip))) stop("item parameters not a matrix")
 if (nrow(ip) != ncol(resp)) stop("responses - item parameters mismatch")
 np = nrow(resp)
 o = sapply(1:np, function(i) dpv.one(resp=resp[i,], ip=ip, mu=mu, s=sigma, n=n))
 return(t(o)) 
}
