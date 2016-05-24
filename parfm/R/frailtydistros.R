
################################################################################
#  Frailty distributions                                                       #
################################################################################
#                                                                              #
#  These are functions with parameters                                         #
#   - k        : the order of the derivative of the Laplace transform          #
#   - s        : the argument of the Laplace transform                         #
#   - theta/nu/sigma2
#              : the heterogeneity parameter of the frailty distribution       #
#   - what     : the quantity to be returned by the function,                  #
#                either "logLT" for \log[ (-1)^k \mathcal L^(k)(s) ]           # 
#                with \mathcal L(s) the Laplace transofrm                      #
#                and \mathcal L^(k)(s) its k-th derivative,                    #
#                or "tau", the Kendall's Tau                                   #
#                                                                              #
#   - correct  : (only for possta) the correction to use in case of many       #
#                events per cluster to get finite likelihood values.           #
#                When correct!=0 the likelihood is divided by                  #
#                10^(#clusters * correct) for computation,                     # 
#                but the value of the log-likelihood in the output             #
#                is the re-adjusted value.                                     #
#                                                                              #
#   Date: December, 19, 2011                                                   #
#   Last modification on: May 1, 2013                                          #
################################################################################



################################################################################
#                                                                              #
#   No frailty distribution                                                    #
#                                                                              #
#                                                                              #
#   Date: December 21, 2011                                                    #
#   Last modification on: December 27, 2011                                    #
################################################################################

fr.none <- function(s,
                    what="logLT"){
  if (what=="logLT")
    return(-s)
  else if (what == "tau")
    return(NULL)
}



################################################################################
#                                                                              #
#   Gamma frailty distribution                                                 #
#                                                                              #
#   Density:                                                                   #
#    f(u) = \frac{                                                             #
#     \theta^{-\frac1\theta}  u^{\frac1\theta - 1}  \exp( -u / \theta)         #
#    }{                                                                        #
#     \Gamma(1 / \theta) }                                                     #
#                                                                              #
#   Arguments of fr.gamma:                                                     #
#     [1] k = 0, 1, ...                                                        #
#     [2] s > 0                                                                #
#     [3] theta > 0                                                            #
#                                                                              #
#   Date: December 21, 2011                                                    #
#   Last modification on: December 27, 2011                                    #
################################################################################

fr.gamma <- function(k,
                     s, 
                     theta, 
                     what="logLT"){
  if (what=="logLT") {
    res <- ifelse(k == 0, 
                  - 1 / theta  * log(1 + theta * s),
                  - (k + 1 / theta) * log(1 + theta * s) +
                    sum(log(1 + (seq(from=0, to=k-1, by=1) * theta))))
    return(res)
  }
  else if (what == "tau")
    return(theta / (theta + 2))
}



################################################################################
#                                                                              #
#   Inverse Gaussian frailty distribution                                      #
#                                                                              #
#   Density:                                                                   #
#    f(u) = \frac1{ \sqrt{2 \pi \theta} }  u^{-\frac32}                        #
#           \exp( -\frac{(u-1)^2}{2 \theta u} )                                #
#                                                                              #
#   Arguments of fr.ingau:                                                     #
#     [1] k = 0, 1, ...                                                        #
#     [2] s > 0                                                                #
#     [3] theta > 0                                                            #
#                                                                              #
#   Date: December 20, 2011                                                    #
#   Last modification on: January 13, 2012                                     #
################################################################################

fr.ingau <- function(k, 
                     s, 
                     theta, 
                     what="logLT"){
  if (what=="logLT") {
    # separate sqrt's for numerical reasons in case of very small theta!
    z <- theta^(-0.5) * sqrt(2 * s + theta^(-1))
    res <- ifelse(k == 0,
                  1 / theta * (1 - sqrt(1 + 2 * theta * s)),
                  - k / 2 * log(2 * theta * s + 1) +
                    log(besselK(z, k - 0.5)) - 
                    (log(pi / (2 * z)) / 2 - z) +
                    1 / theta * (1 - sqrt(1 + 2 * theta * s)))
    return(res)
  }
  else if (what == "tau") {
    integrand <- function(u) {
      return(exp(-u) / u)
    }
    int <- integrate(integrand,
                     lower=(2/theta), 
                     upper=Inf)$value
    tau <- 0.5 - (1 / theta) + (2 * theta^(-2) * exp(2 / theta) * int)
    if (is.nan(tau) || tau < 0)
      tau <- paste("The value of 'theta' is too small",
                   "for computing the Kendall's Tau numerically!")
    return(tau)
  }
}



################################################################################
#                                                                              #
#   Sum of the polynomials Omega for the Positive Stable frailty distribution  #
#                                                                              #
#                                                                              #
#   Date: December 20, 2011                                                    #
#   Last modification on: January 10, 2012                                     #
################################################################################

J <- function(k, s, nu, Omega, correct){  
  if(k == 0) sum <- 10^-correct else {
    sum <- 0
    for(m in 0:(k - 1)) {
      sum <- sum + (Omega[k, m + 1] * s^(-m * (1 - nu)))
    }
  }
  return(sum)
}



################################################################################
#                                                                              #
#   Positive Stable frailty distribution                                       #
#                                                                              #
#   Density:                                                                   #
#    f(u) = -\frac1{\pi u}                                                     #
#           \sum_{k=1}^\infty \frac{ \Gamma( k (1 - \nu) + 1) }{ k! }          #
#           ( -u^{\nu-1} )^k  \sin( (1-\nu) k \pi)                             #
#                                                                              #
#   Arguments of fr.possta:                                                    #
#     [1] k = 0, 1, ...                                                        #
#     [2] s > 0                                                                #
#     [3] nu in (0, 1)                                                         #
#     [4] Omega is the matrix that contains the omega's                        #
#                                                                              #
#   Date: December 20, 2011                                                    #
#   Last modification on: January 10, 2012                                     #
################################################################################

fr.possta <- function(k,
                      s,
                      nu,
                      Omega,
                      what="logLT",
                      correct){
  if (what=="logLT") {
    res <- k * (log(1 - nu) - nu * log(s)) - s^(1 - nu) + 
      log(J(k, s, nu, Omega, correct))
    return(res)
  }
  else if (what == "tau")
    return(nu)
}



################################################################################
#                                                                              #
#   Lognormal frailty distribution                                             #
#                                                                              #
#   Density:                                                                   #
#    f(u) = -\frac1{\sqrt{2\pi \theta}}                                        #
#           \exp\{ -\frac{u^2}{2\theta}\}                                      #
#                                                                              #
#   Arguments of fr.lognormal:                                                 #
#     [1] k = 0, 1, ...                                                        #
#     [2] s > 0                                                                #
#     [3] sigma > 0                                                            #
#                                                                              #
#   Date: October 16, 2012                                                     #
#   Last modification on: September 2, 2015                                    #
################################################################################
g <- function(w, k, s, sigma2) {
    - k * w + exp(w) * s + .5 * w^2 / sigma2
}

g1 <- function(w, k, s, sigma2) {
    - k + exp(w) * s + w / sigma2
}

g2 <- function(w, k, s, sigma2) {
    exp(w) * s + 1 / sigma2    
} 

Lapl <- Vectorize(function(s, k, sigma2) {
    # Find wTilde = max(g(w)) so that g'(wTilde; k, s, theta) = 0
    WARN <- getOption("warn")
    options(warn=-1)
    wTilde <- optimize(f=g, c(-1e10, 1e10), maximum=FALSE,
                       k=k, s=s, sigma2=sigma2)$minimum
    options(warn=WARN)
    
    # Approximate the integral via Laplacian method
    res <- (-1)^k * 
        exp(-g(w=wTilde, k=k, s=s, sigma2=sigma2)) /
        sqrt(sigma2 * g2(w=wTilde, k=k, s=s, sigma2=sigma2))
    return(res)
}, 's')

fr.lognormal <- function(k,
                         s,
                         sigma2,
                         what="logLT"){
    if (what=="logLT") {
        # Find wTilde = max(g(w)) so that g'(wTilde; k, s, theta) = 0
        WARN <- getOption("warn")
        options(warn=-1)
        wTilde <- optimize(f=g, c(-1e10, 1e10), maximum=FALSE,
                           k=k, s=s, sigma2=sigma2)$minimum
        options(warn=WARN)
        
        # Approximate the integral via Laplacian method
        res <- (-1)^k * 
            exp(-g(w=wTilde, k=k, s=s, sigma2=sigma2)) /
            sqrt(sigma2 * g2(w=wTilde, k=k, s=s, sigma2=sigma2))
        return(res)
    }
    else if (what == "tau") {
        intTau <- Vectorize(function(x, intTau.sigma2=sigma2) {
            res <- x * 
                Lapl(s=x, k = 0, sigma2 = intTau.sigma2) *
                Lapl(s=x, k = 2, sigma2 = intTau.sigma2)
            return(res)
        }, "x")
        
        tauRes <- 4 * integrate(f=intTau, lower=0, upper=Inf, intTau.sigma2=sigma2)$value - 1
        return(tauRes)
    }
}