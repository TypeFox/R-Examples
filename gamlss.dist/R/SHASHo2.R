###############################
##### SHASHo2 DISTRIBUTION #####
###############################

# This is the original SHASH distribution from Jones and Pewse3y (2009) Biometrika. 1-20
# see page 2 equation (2) of the paper (divided by sigma for the density of the unstandardized variable)
# reparameterized to a new sigma = old sigma/tau  as recommended for fitting on page 8 of the paper
#--------------------------------------------------------------------------------
##Gamlss family function:
SHASHo2 <- function (mu.link = "identity", sigma.link = "log", 
                    nu.link = "identity", tau.link = "log") {
    mstats <- checklink("mu.link", "Sinh-Arcsinh", substitute(mu.link), 
        c("inverse", "log", "identity", "own"))
    dstats <- checklink("sigma.link", "Sinh-Arcsinh", substitute(sigma.link), 
        c("inverse", "log", "identity", "own"))
    vstats <- checklink("nu.link", "Sinh-Arcsinh", substitute(nu.link), 
        c("inverse", "log", "identity", "own"))
    tstats <- checklink("tau.link", "Sinh-Arcsinh", substitute(tau.link), 
        c("inverse", "log", "identity", "own"))
    structure(list(family = c("SHASHo2", "Sinh-Arcsinh"), 
        parameters = list(mu = TRUE, sigma = TRUE, nu = TRUE, tau = TRUE), 
        nopar = 4, type = "Continuous", 
        mu.link = as.character(substitute(mu.link)), 
        sigma.link = as.character(substitute(sigma.link)), 
        nu.link = as.character(substitute(nu.link)), 
        tau.link = as.character(substitute(tau.link)), 
        mu.linkfun = mstats$linkfun, sigma.linkfun = dstats$linkfun, 
        nu.linkfun = vstats$linkfun, tau.linkfun = tstats$linkfun, 
        mu.linkinv = mstats$linkinv, sigma.linkinv = dstats$linkinv, 
        nu.linkinv = vstats$linkinv, tau.linkinv = tstats$linkinv, 
        mu.dr = mstats$mu.eta, sigma.dr = dstats$mu.eta,      
        nu.dr = vstats$mu.eta, tau.dr = tstats$mu.eta,                    
        dldm = function(y, mu, sigma, nu, tau) {
    z <- (y-mu)/(sigma*tau)
   c <- cosh(tau*asinh(z)-nu)
   r <- sinh(tau*asinh(z)-nu)
          x <- tau*asinh(z)-nu
          dldm <- (1/(tau*sigma*(1+z^2)^(1/2)))*(r*tau*c - ((tau*sinh(x))/c)  +
                   z/(1+z^2)^(1/2))
          dldm
        }, 
        d2ldm2 = function(y, mu, sigma, nu, tau) {
   z <- (y-mu)/(sigma*tau)
   c <- cosh(tau*asinh(z)-nu)
   r <- sinh(tau*asinh(z)-nu)
          x <- tau*asinh(z)-nu
          dldm <- (1/(tau*sigma*(1+z^2)^(1/2)))*(r*tau*c - ((tau*sinh(x))/c)  +
                   z/(1+z^2)^(1/2))
          d2ldm2 <- -dldm*dldm
          d2ldm2
        }, 
        dldd = function(y, mu, sigma, nu, tau) {
        #   dd <- numeric.deriv(dSHASHo(y, mu, sigma, 
       #         nu, tau, log = TRUE), "sigma", delta=0.001)
       #   dldd <- as.vector(attr(dd, "gradient"))
          
   z <- (y-mu)/(sigma*tau)
   c <- cosh(tau*asinh(z)-nu)
   r <- sinh(tau*asinh(z)-nu)
          x <- tau*asinh(z)-nu
          dldd <- ((z)/(sigma*(1+z^2)^(1/2)))*(r*tau*c-((tau*sinh(x))/c)+
                  z/((1+z^2)^(1/2)))-(1/sigma)
          dldd
        }, 
        d2ldd2 = function(y, mu, sigma, nu, tau) {
  #        dd <- numeric.deriv(dSHASHo2(y, mu, sigma, 
  #              nu, tau, log = TRUE), "sigma", delta=0.001)
  #        dldd <- as.vector(attr(dd, "gradient"))
  
  z <- (y-mu)/(sigma*tau)
   c <- cosh(tau*asinh(z)-nu)
   r <- sinh(tau*asinh(z)-nu)
          x <- tau*asinh(z)-nu
          dldd <- ((z)/(sigma*(1+z^2)^(1/2)))*(r*tau*c-((tau*sinh(x))/c)+
                  z/((1+z^2)^(1/2)))-(1/sigma)
          d2ldd2 <- -dldd*dldd
          d2ldd2
        }, 
        dldv = function(y, mu, sigma, nu, tau) {
 z <- (y-mu)/(sigma*tau)
   c <- cosh(tau*asinh(z)-nu)
   r <- sinh(tau*asinh(z)-nu)
            dldv<- r*cosh(tau*asinh(z)-nu)-(1/c)*sinh(tau*asinh(z)-nu)
 #         nd <- numeric.deriv(dSHASHo2(y, mu, sigma, 
 #               nu, tau, log = TRUE), "nu", delta=0.001)
 #         dldv <- as.vector(attr(nd, "gradient"))
          dldv
        }, 
        d2ldv2 = function(y, mu, sigma, nu, tau) {
z <- (y-mu)/(sigma*tau)
   c <- cosh(tau*asinh(z)-nu)
   r <- sinh(tau*asinh(z)-nu)
            dldv<- r*cosh(tau*asinh(z)-nu)-(1/c)*sinh(tau*asinh(z)-nu)
 #         nd <- numeric.deriv(dSHASHo2(y, mu, sigma, 
 #               nu, tau, log = TRUE), "nu", delta=0.001)
 #         dldv <- as.vector(attr(nd, "gradient"))
          d2ldv2 <- -dldv*dldv
          d2ldv2
        }, 
        dldt = function(y, mu, sigma, nu, tau) {
z <- (y-mu)/(sigma*tau)
   c <- cosh(tau*asinh(z)-nu)
   r <- sinh(tau*asinh(z)-nu)
          dldt <- (-r*cosh(tau*asinh(z)-nu)+(1/c)*sinh(tau*asinh(z)-nu))*
                  (asinh(z)+(1/(1+z^2)^(1/2))*(-z)) + (z^2/(tau*(1+z^2)))
         
 #                 td <- numeric.deriv(dSHASHo2(y, mu, sigma, 
 #               nu, tau, log = TRUE), "tau", delta=0.001)
 #         dldt <- as.vector(attr(td, "gradient"))
 
          dldt
        }, 
        d2ldt2 = function(y, mu, sigma, nu, tau) {
z <- (y-mu)/(sigma*tau)
   c <- cosh(tau*asinh(z)-nu)
   r <- sinh(tau*asinh(z)-nu)
          dldt <- (-r*cosh(tau*asinh(z)-nu)+(1/c)*sinh(tau*asinh(z)-nu))*
                  (asinh(z)+(1/(1+z^2)^(1/2))*(-z)) + (z^2/(tau*(1+z^2)))
         
 #                 td <- numeric.deriv(dSHASHo2(y, mu, sigma, 
 #               nu, tau, log = TRUE), "tau", delta=0.001)
 #         dldt <- as.vector(attr(td, "gradient"))
          d2ldt2 <- -dldt*dldt
          d2ldt2
        }, 
        d2ldmdd = function(y, mu, sigma, nu, tau) {
     z <- (y-mu)/(sigma*tau)
   c <- cosh(tau*asinh(z)-nu)
   r <- sinh(tau*asinh(z)-nu)
          x <- tau*asinh(z)-nu
          dldm <- (1/(tau*sigma*(1+z^2)^(1/2)))*(r*tau*c - ((tau*sinh(x))/c)  +
                   z/(1+z^2)^(1/2))
  
          dldd <- ((z)/(sigma*(1+z^2)^(1/2)))*(r*tau*c-((tau*sinh(x))/c)+
                  z/((1+z^2)^(1/2)))-(1/sigma)
          d2ldmdd <- -dldd*dldm
          d2ldmdd
        }, 
        d2ldmdv = function(y, mu, sigma, nu, tau) {
    z <- (y-mu)/(sigma*tau)
   c <- cosh(tau*asinh(z)-nu)
   r <- sinh(tau*asinh(z)-nu)
          x <- tau*asinh(z)-nu
          dldm <- (1/(tau*sigma*(1+z^2)^(1/2)))*(r*tau*c - ((tau*sinh(x))/c)  +
                   z/(1+z^2)^(1/2))
           dldv<- r*cosh(tau*asinh(z)-nu)-(1/c)*sinh(tau*asinh(z)-nu)
          d2ldmdv <- -dldm*dldv
          d2ldmdv
        }, 
        d2ldmdt = function(y, mu, sigma, nu, tau) {
    z <- (y-mu)/(sigma*tau)
   c <- cosh(tau*asinh(z)-nu)
   r <- sinh(tau*asinh(z)-nu)
          x <- tau*asinh(z)-nu
          dldm <- (1/(tau*sigma*(1+z^2)^(1/2)))*(r*tau*c - ((tau*sinh(x))/c)  +
                   z/(1+z^2)^(1/2))
        dldt <- (-r*cosh(tau*asinh(z)-nu)+(1/c)*sinh(tau*asinh(z)-nu))*
                  (asinh(z)+(1/(1+z^2)^(1/2))*(-z)) + (z^2/(tau*(1+z^2)))
          d2ldmdt <- -dldm*dldt
          d2ldmdt
        }, 
        d2ldddv = function(y, mu, sigma, nu, tau) {
   z <- (y-mu)/(sigma*tau)
   c <- cosh(tau*asinh(z)-nu)
   r <- sinh(tau*asinh(z)-nu)
          x <- tau*asinh(z)-nu
          dldd <- ((z)/(sigma*(1+z^2)^(1/2)))*(r*tau*c-((tau*sinh(x))/c)+
                  z/((1+z^2)^(1/2)))-(1/sigma)
           dldv<- r*cosh(tau*asinh(z)-nu)-(1/c)*sinh(tau*asinh(z)-nu)
          d2ldddv <- -dldd*dldv
          d2ldddv
        }, 
        d2ldddt = function(y, mu, sigma, nu, tau) {
   z <- (y-mu)/(sigma*tau)
   c <- cosh(tau*asinh(z)-nu)
   r <- sinh(tau*asinh(z)-nu)
          x <- tau*asinh(z)-nu
          dldd <- ((z)/(sigma*(1+z^2)^(1/2)))*(r*tau*c-((tau*sinh(x))/c)+
                  z/((1+z^2)^(1/2)))-(1/sigma)
          
          dldt <- (-r*cosh(tau*asinh(z)-nu)+(1/c)*sinh(tau*asinh(z)-nu))*
                  (asinh(z)+(1/(1+z^2)^(1/2))*(-z)) + (z^2/(tau*(1+z^2)))
          d2ldddt <- -dldd*dldt
          d2ldddt
        }, 
        d2ldvdt = function(y, mu, sigma, nu, tau) {
        z <- (y-mu)/(sigma*tau)
   c <- cosh(tau*asinh(z)-nu)
   r <- sinh(tau*asinh(z)-nu)
                  dldv<- r*cosh(tau*asinh(z)-nu)-(1/c)*sinh(tau*asinh(z)-nu)
         dldt <- (-r*cosh(tau*asinh(z)-nu)+(1/c)*sinh(tau*asinh(z)-nu))*
                  (asinh(z)+(1/(1+z^2)^(1/2))*(-z)) + (z^2/(tau*(1+z^2)))
          d2ldvdt <- -dldv*dldt
          d2ldvdt 
        }, 
        G.dev.incr = function(y, mu, sigma, nu, tau, ...) -2 * 
          dSHASHo2(y, mu, sigma, nu, tau, log = TRUE), 
        rqres = expression(rqres(pfun = "pSHASHo2", 
        type = "Continuous", y = y, mu = mu, sigma = sigma, nu = nu, tau = tau)), 
        mu.initial = expression(mu <- (y + mean(y))/2), 
        sigma.initial = expression(sigma <- rep(sd(y)/5, length(y))), 
        nu.initial = expression(nu <- rep(0.5, length(y))), 
        tau.initial = expression(tau <- rep(0.5, length(y))), 
        mu.valid = function(mu) TRUE, 
        sigma.valid = function(sigma) all(sigma > 0), 
        nu.valid = function(nu) TRUE, 
        tau.valid = function(tau) all(tau > 0), 
        y.valid = function(y) TRUE), 
        class = c("gamlss.family", "family"))
}
#-------------------------------------------------------------------------------
#Probability density function
dSHASHo2 <- function(x, mu=0, sigma=1, nu=0, tau=1, log = FALSE){
  if (any(sigma < 0)) 
      stop(paste("sigma must be positive", "\n", ""))
  if (any(tau < 0)) 
      stop(paste("tau must be positive", "\n", ""))
#  sigmat<- sigma/tau
#  z <- (x-mu)/(tau*sigmat)
#  c <- (1/2)*(exp(tau*asinh(z)-nu)+exp(-(tau*asinh(z)-nu)))
#  r <- (1/2) * (exp(tau * asinh(z) - nu) - exp(-tau * asinh(z) + nu))
#  loglik <- -log(sigma) -(1/2)*log(2*pi)-(1/2)*log(1+z^2)+
#            log(tau) + log(c) - (1/2)*r^2
   z <- (x-mu)/(sigma*tau)
   c <- cosh(tau*asinh(z)-nu)
   r <- sinh(tau*asinh(z)-nu)
   loglik <- -log(sigma) -0.5*log(2*pi) -0.5*log(1+(z^2)) +log(c) -0.5*(r^2)
  
  if (log == FALSE)
    fy <- exp(loglik)
  else fy <- loglik
  fy
}
#-------------------------------------------------------------------------------
#Cumulative density function
pSHASHo2 <- function(q, mu=0, sigma=1, nu=0, tau=1, 
                    lower.tail = TRUE, log.p = FALSE){
  if (any(sigma < 0)) 
      stop(paste("sigma must be positive", "\n", ""))
  if (any(tau < 0)) 
      stop(paste("tau must be positive", "\n", ""))
 
 z <- (q-mu)/(sigma*tau)
   c <- cosh(tau*asinh(z)-nu)
   r <- sinh(tau*asinh(z)-nu)
  p <- pNO(r)
  if (lower.tail == TRUE) 
        p <- p
  else p <- 1 - p
  if (log.p == FALSE) 
        p <- p
  else p <- log(p)
  p
}
#-------------------------------------------------------------------------------
#Quantile function
qSHASHo2 <- function(p, mu=0, sigma=1, nu=0, tau=1, lower.tail = TRUE, 
                    log.p = FALSE){
   if (log.p==TRUE) p <- exp(p) else p <- p
   if (any(p <= 0)|any(p >= 1))  stop(paste("p must be between 0 and 1", "\n", ""))       
   if (lower.tail==TRUE) p <- p else p <- 1-p
   
   y <- mu + (tau*sigma)*sinh((1/tau)*asinh(qnorm(p))+(nu/tau))
   y
}
#--------------------------------------------------------------------------------
#Random generation
rSHASHo2 <- function(n, mu=0, sigma=1, nu=0, tau=1){
   if (any(n <= 0)) 
        stop(paste("n must be a positive integer", "\n", ""))  
   n <- ceiling(n)
   p <- runif(n)
   r <- qSHASHo2(p, mu = mu, sigma = sigma, nu=nu, tau=tau)
   r 
}
