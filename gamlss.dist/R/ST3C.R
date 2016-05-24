# updated 27_11_2007
ST3C <- function (mu.link="identity", sigma.link="log", nu.link ="log", tau.link="log")
{
    mstats <- checklink("mu.link", "skew t type 2", substitute(mu.link), 
                         c("inverse", "log", "identity", "own"))
    dstats <- checklink("sigma.link", "skew t type 2", substitute(sigma.link), 
                         c("inverse", "log", "identity", "own"))
    vstats <- checklink("nu.link", "skew t type 2",substitute(nu.link), 
                         c("inverse", "log", "identity", "own"))
    tstats <- checklink("tau.link", "skew t type 2 ",substitute(tau.link), 
                         c("inverse", "log", "identity", "own")) 
    structure(
          list(family = c("ST3",  "skew t type 3"),
           parameters = list(mu=TRUE, sigma=TRUE, nu=TRUE, tau=TRUE), 
                nopar = 4, 
                 type = "Continuous",
              mu.link = as.character(substitute(mu.link)),  
           sigma.link = as.character(substitute(sigma.link)), 
              nu.link = as.character(substitute(nu.link)), 
             tau.link = as.character(substitute(tau.link)), 
           mu.linkfun = mstats$linkfun, 
        sigma.linkfun = dstats$linkfun, 
           nu.linkfun = vstats$linkfun,
           tau.linkfun = tstats$linkfun,  
           mu.linkinv = mstats$linkinv, 
        sigma.linkinv = dstats$linkinv,
           nu.linkinv = vstats$linkinv,
           tau.linkinv = tstats$linkinv, 
                mu.dr = mstats$mu.eta, 
             sigma.dr = dstats$mu.eta, 
                nu.dr = vstats$mu.eta,
               tau.dr = tstats$mu.eta, 
                 dldm = function(y,mu,sigma,nu,tau){
					 n = c(length(y), length(mu), length(sigma), length(nu), length(tau))
					 maxn = max(n)
					 if(n[1]!=maxn) y = rep(y[1], maxn)
					 if(n[2]!=maxn) mu    = rep(mu[1], maxn)
					 if(n[3]!=maxn) sigma = rep(sigma[1], maxn)
					 if(n[4]!=maxn) nu  = rep(nu[1], maxn)
					 if(n[5]!=maxn) tau = rep(tau[1], maxn)
					 ans = double(maxn)
					 sol = try(.C("c_st3_dldm", y = as.double(y), mu = as.double(mu), 
									 sigma = as.double(sigma), nu = as.double(nu), 
									 tau = as.double(tau), ans = ans, n = as.integer(maxn), 
									 PACKAGE="gamlss.dist"), silent = TRUE)
					 if(inherits(sol, 'try-error')){
						 return(sol)
					 } else{
						 return(sol$ans)
					 }                        
				 },
       		 	d2ldm2 = function(y,mu,sigma,nu,tau){
					n = c(length(y), length(mu), length(sigma), length(nu), length(tau))
					maxn = max(n)
					if(n[1]!=maxn) y = rep(y[1], maxn)
					if(n[2]!=maxn) mu    = rep(mu[1], maxn)
					if(n[3]!=maxn) sigma = rep(sigma[1], maxn)
					if(n[4]!=maxn) nu  = rep(nu[1], maxn)
					if(n[5]!=maxn) tau = rep(tau[1], maxn)
					ans = double(maxn)
					sol = try(.C("c_st3_d2ldm2", y = as.double(y), mu = as.double(mu), 
									sigma = as.double(sigma), nu = as.double(nu), 
									tau = as.double(tau), ans = ans, n = as.integer(maxn), 
									PACKAGE="gamlss.dist"), silent = TRUE)
					if(inherits(sol, 'try-error')){
						return(sol)
					} else{
						return(sol$ans)
					}                        
				},
                 dldd = function(y,mu,sigma,nu,tau) {
					 n = c(length(y), length(mu), length(sigma), length(nu), length(tau))
					 maxn = max(n)
					 if(n[1]!=maxn) y = rep(y[1], maxn)
					 if(n[2]!=maxn) mu    = rep(mu[1], maxn)
					 if(n[3]!=maxn) sigma = rep(sigma[1], maxn)
					 if(n[4]!=maxn) nu  = rep(nu[1], maxn)
					 if(n[5]!=maxn) tau = rep(tau[1], maxn)
					 ans = double(maxn)
					 sol = try(.C("c_st3_dldd", y = as.double(y), mu = as.double(mu), 
									 sigma = as.double(sigma), nu = as.double(nu), 
									 tau = as.double(tau), ans = ans, n = as.integer(maxn), 
									 PACKAGE="gamlss.dist"), silent = TRUE)
					 if(inherits(sol, 'try-error')){
						 return(sol)
					 } else{
						 return(sol$ans)
					 }                        
				 },
				 d2ldd2 = function(y,mu,sigma,nu,tau) {
					 n = c(length(y), length(mu), length(sigma), length(nu), length(tau))
					 maxn = max(n)
					 if(n[1]!=maxn) y = rep(y[1], maxn)
					 if(n[2]!=maxn) mu    = rep(mu[1], maxn)
					 if(n[3]!=maxn) sigma = rep(sigma[1], maxn)
					 if(n[4]!=maxn) nu  = rep(nu[1], maxn)
					 if(n[5]!=maxn) tau = rep(tau[1], maxn)
					 ans = double(maxn)
					 sol = try(.C("c_st3_d2ldd2", y = as.double(y), mu = as.double(mu), 
									 sigma = as.double(sigma), nu = as.double(nu), 
									 tau = as.double(tau), ans = ans, n = as.integer(maxn), 
									 PACKAGE="gamlss.dist"), silent = TRUE)
					 if(inherits(sol, 'try-error')){
						 return(sol)
					 } else{
						 return(sol$ans)
					 }                        
				 },
                 dldv = function(y,mu,sigma,nu,tau) {
					 n = c(length(y), length(mu), length(sigma), length(nu), length(tau))
					 maxn = max(n)
					 if(n[1]!=maxn) y = rep(y[1], maxn)
					 if(n[2]!=maxn) mu    = rep(mu[1], maxn)
					 if(n[3]!=maxn) sigma = rep(sigma[1], maxn)
					 if(n[4]!=maxn) nu  = rep(nu[1], maxn)
					 if(n[5]!=maxn) tau = rep(tau[1], maxn)
					 ans = double(maxn)
					 sol = try(.C("c_st3_dldv", y = as.double(y), mu = as.double(mu), 
									 sigma = as.double(sigma), nu = as.double(nu), 
									 tau = as.double(tau), ans = ans, n = as.integer(maxn), 
									 PACKAGE="gamlss.dist"), silent = TRUE)
					 if(inherits(sol, 'try-error')){
						 return(sol)
					 } else{
						 return(sol$ans)
					 }                
				 },
               d2ldv2 = function(y,mu,sigma,nu,tau) { 
				   n = c(length(y), length(mu), length(sigma), length(nu), length(tau))
				   maxn = max(n)
				   if(n[1]!=maxn) y = rep(y[1], maxn)
				   if(n[2]!=maxn) mu    = rep(mu[1], maxn)
				   if(n[3]!=maxn) sigma = rep(sigma[1], maxn)
				   if(n[4]!=maxn) nu  = rep(nu[1], maxn)
				   if(n[5]!=maxn) tau = rep(tau[1], maxn)
				   ans = double(maxn)
				   sol = try(.C("c_st3_d2ldv2", y = as.double(y), mu = as.double(mu), 
								   sigma = as.double(sigma), nu = as.double(nu), 
								   tau = as.double(tau), ans = ans, n = as.integer(maxn), 
								   PACKAGE="gamlss.dist"), silent = TRUE)
				   if(inherits(sol, 'try-error')){
					   return(sol)
				   } else{
					   return(sol$ans)
				   }                
			   },
                 dldt = function(y,mu,sigma,nu,tau) { 
					 n = c(length(y), length(mu), length(sigma), length(nu), length(tau))
					 maxn = max(n)
					 if(n[1]!=maxn) y = rep(y[1], maxn)
					 if(n[2]!=maxn) mu    = rep(mu[1], maxn)
					 if(n[3]!=maxn) sigma = rep(sigma[1], maxn)
					 if(n[4]!=maxn) nu  = rep(nu[1], maxn)
					 if(n[5]!=maxn) tau = rep(tau[1], maxn)
					 ans = double(maxn)
					 sol = try(.C("c_st3_dldt", y = as.double(y), mu = as.double(mu), 
									 sigma = as.double(sigma), nu = as.double(nu), 
									 tau = as.double(tau), ans = ans, n = as.integer(maxn), 
									 PACKAGE="gamlss.dist"), silent = TRUE)
					 if(inherits(sol, 'try-error')){
						 return(sol)
					 } else{
						 return(sol$ans)
					 }                
				 },
               d2ldt2 = function(y,mu,sigma,nu,tau) {
				   n = c(length(y), length(mu), length(sigma), length(nu), length(tau))
				   maxn = max(n)
				   if(n[1]!=maxn) y = rep(y[1], maxn)
				   if(n[2]!=maxn) mu    = rep(mu[1], maxn)
				   if(n[3]!=maxn) sigma = rep(sigma[1], maxn)
				   if(n[4]!=maxn) nu  = rep(nu[1], maxn)
				   if(n[5]!=maxn) tau = rep(tau[1], maxn)
				   ans = double(maxn)
				   sol = try(.C("c_st3_d2ldt2", y = as.double(y), mu = as.double(mu), 
								   sigma = as.double(sigma), nu = as.double(nu), 
								   tau = as.double(tau), ans = ans, n = as.integer(maxn), 
								   PACKAGE="gamlss.dist"), silent = TRUE)
				   if(inherits(sol, 'try-error')){
					   return(sol)
				   } else{
					   return(sol$ans)
				   }                
			   },
              	d2ldmdd = function(y,mu,sigma,nu,tau){
					n = c(length(y), length(mu), length(sigma), length(nu), length(tau))
					maxn = max(n)
					if(n[1]!=maxn) y = rep(y[1], maxn)
					if(n[2]!=maxn) mu    = rep(mu[1], maxn)
					if(n[3]!=maxn) sigma = rep(sigma[1], maxn)
					if(n[4]!=maxn) nu  = rep(nu[1], maxn)
					if(n[5]!=maxn) tau = rep(tau[1], maxn)
					ans = double(maxn)
					sol = try(.C("c_st3_d2ldmdd", y = as.double(y), mu = as.double(mu), 
									sigma = as.double(sigma), nu = as.double(nu), 
									tau = as.double(tau), ans = ans, n = as.integer(maxn), 
									PACKAGE="gamlss.dist"), silent = TRUE)
					if(inherits(sol, 'try-error')){
						return(sol)
					} else{
						return(sol$ans)
					}                
				},
              d2ldmdv = function(y,mu,sigma,nu,tau){
				  n = c(length(y), length(mu), length(sigma), length(nu), length(tau))
				  maxn = max(n)
				  if(n[1]!=maxn) y = rep(y[1], maxn)
				  if(n[2]!=maxn) mu    = rep(mu[1], maxn)
				  if(n[3]!=maxn) sigma = rep(sigma[1], maxn)
				  if(n[4]!=maxn) nu  = rep(nu[1], maxn)
				  if(n[5]!=maxn) tau = rep(tau[1], maxn)
				  ans = double(maxn)
				  sol = try(.C("c_st3_d2ldmdv", y = as.double(y), mu = as.double(mu), 
								  sigma = as.double(sigma), nu = as.double(nu), 
								  tau = as.double(tau), ans = ans, n = as.integer(maxn), 
								  PACKAGE="gamlss.dist"), silent = TRUE)
				  if(inherits(sol, 'try-error')){
					  return(sol)
				  } else{
					  return(sol$ans)
				  }                
			  },
              d2ldmdt = function(y,mu,sigma,nu,tau){  
				  n = c(length(y), length(mu), length(sigma), length(nu), length(tau))
				  maxn = max(n)
				  if(n[1]!=maxn) y = rep(y[1], maxn)
				  if(n[2]!=maxn) mu    = rep(mu[1], maxn)
				  if(n[3]!=maxn) sigma = rep(sigma[1], maxn)
				  if(n[4]!=maxn) nu  = rep(nu[1], maxn)
				  if(n[5]!=maxn) tau = rep(tau[1], maxn)
				  ans = double(maxn)
				  sol = try(.C("c_st3_d2ldmdt", y = as.double(y), mu = as.double(mu), 
								  sigma = as.double(sigma), nu = as.double(nu), 
								  tau = as.double(tau), ans = ans, n = as.integer(maxn), 
								  PACKAGE="gamlss.dist"), silent = TRUE)
				  if(inherits(sol, 'try-error')){
					  return(sol)
				  } else{
					  return(sol$ans)
				  }                
			  },
              d2ldddv = function(y,mu,sigma,nu,tau){
				  n = c(length(y), length(mu), length(sigma), length(nu), length(tau))
				  maxn = max(n)
				  if(n[1]!=maxn) y = rep(y[1], maxn)
				  if(n[2]!=maxn) mu    = rep(mu[1], maxn)
				  if(n[3]!=maxn) sigma = rep(sigma[1], maxn)
				  if(n[4]!=maxn) nu  = rep(nu[1], maxn)
				  if(n[5]!=maxn) tau = rep(tau[1], maxn)
				  ans = double(maxn)
				  sol = try(.C("c_st3_d2ldddv", y = as.double(y), mu = as.double(mu), 
								  sigma = as.double(sigma), nu = as.double(nu), 
								  tau = as.double(tau), ans = ans, n = as.integer(maxn), 
								  PACKAGE="gamlss.dist"), silent = TRUE)
				  if(inherits(sol, 'try-error')){
					  return(sol)
				  } else{
					  return(sol$ans)
				  }                
			  },
              d2ldddt = function(y,mu,sigma,nu,tau){
				  n = c(length(y), length(mu), length(sigma), length(nu), length(tau))
				  maxn = max(n)
				  if(n[1]!=maxn) y = rep(y[1], maxn)
				  if(n[2]!=maxn) mu    = rep(mu[1], maxn)
				  if(n[3]!=maxn) sigma = rep(sigma[1], maxn)
				  if(n[4]!=maxn) nu  = rep(nu[1], maxn)
				  if(n[5]!=maxn) tau = rep(tau[1], maxn)
				  ans = double(maxn)
				  sol = try(.C("c_st3_d2ldddt", y = as.double(y), mu = as.double(mu), 
								  sigma = as.double(sigma), nu = as.double(nu), 
								  tau = as.double(tau), ans = ans, n = as.integer(maxn), 
								  PACKAGE="gamlss.dist"), silent = TRUE)
				  if(inherits(sol, 'try-error')){
					  return(sol)
				  } else{
					  return(sol$ans)
				  }                
			  },
              d2ldvdt = function(y,mu,sigma,nu,tau){
				  n = c(length(y), length(mu), length(sigma), length(nu), length(tau))
				  maxn = max(n)
				  if(n[1]!=maxn) y = rep(y[1], maxn)
				  if(n[2]!=maxn) mu    = rep(mu[1], maxn)
				  if(n[3]!=maxn) sigma = rep(sigma[1], maxn)
				  if(n[4]!=maxn) nu  = rep(nu[1], maxn)
				  if(n[5]!=maxn) tau = rep(tau[1], maxn)
				  ans = double(maxn)
				  sol = try(.C("c_st3_d2ldvdt", y = as.double(y), mu = as.double(mu), 
								  sigma = as.double(sigma), nu = as.double(nu), 
								  tau = as.double(tau), ans = ans, n = as.integer(maxn), 
								  PACKAGE="gamlss.dist"), silent = TRUE)
				  if(inherits(sol, 'try-error')){
					  return(sol)
				  } else{
					  return(sol$ans)
				  }                
			  },
          G.dev.incr  = function(y,mu,sigma,nu,tau,...) 
                                  -2*dST3(y,mu,sigma,nu,tau,log=TRUE), 
                 rqres = expression(rqres(pfun="pST3", type="Continuous", y=y, mu=mu, sigma=sigma, nu=nu, tau=tau)) ,
            mu.initial = expression(mu <- (y+mean(y))/2), 
         sigma.initial = expression(sigma<- rep(sd(y), length(y))),
            nu.initial = expression(nu <- rep(1, length(y))), 
           tau.initial = expression(tau <-rep(10, length(y))), 
              mu.valid = function(mu) TRUE, 
           sigma.valid = function(sigma)  all(sigma > 0),
              nu.valid = function(nu) all(nu > 0), 
             tau.valid = function(tau) all(tau > 0), 
               y.valid = function(y)  TRUE
          ),
            class = c("gamlss.family","family"))
}
#-----------------------------------------------------------------  
dST3C <- function(x, mu=0, sigma=1, nu=1, tau=10, log=FALSE)
{
      if (any(sigma <= 0))  stop(paste("sigma must be positive", "\n", "")) 
      if (any(nu <= 0))  stop(paste("nu must be positive", "\n", ""))
      if (any(tau <= 0))  stop(paste("tau must be positive", "\n", ""))
	  n = c(length(x), length(mu), length(sigma), length(nu), length(tau))
	  maxn = max(n)
	  if(n[1]!=maxn) x = rep(x[1], maxn)
	  if(n[2]!=maxn) mu    = rep(mu[1], maxn)
	  if(n[3]!=maxn) sigma = rep(sigma[1], maxn)
	  if(n[4]!=maxn) nu  = rep(nu[1], maxn)
	  if(n[5]!=maxn) tau = rep(tau[1], maxn)
	  ans = double(maxn)
	  sol = try(.C("c_st3_dst3", x = as.double(x), mu = as.double(mu), 
					  sigma = as.double(sigma), nu = as.double(nu), 
					  tau = as.double(tau), ans = ans, n = as.integer(maxn), 
					  logr = as.integer(log), PACKAGE="gamlss.dist"), silent = TRUE)
	  if(inherits(sol, 'try-error')){
		  return(sol)
	  } else{
		  return(sol$ans)
	  }
}
#-----------------------------------------------------------------  
pST3C <- function(q, mu=0, sigma=1, nu=1, tau=10, lower.tail = TRUE, log.p = FALSE)
{  
	if (any(sigma <= 0))  stop(paste("sigma must be positive", "\n", "")) 
	if (any(nu <= 0))  stop(paste("nu must be positive", "\n", ""))
	if (any(tau <= 0))  stop(paste("tau must be positive", "\n", ""))
	n = c(length(q), length(mu), length(sigma), length(nu), length(tau))
	maxn = max(n)
	if(n[1]!=maxn) q = rep(q[1], maxn)
	if(n[2]!=maxn) mu    = rep(mu[1], maxn)
	if(n[3]!=maxn) sigma = rep(sigma[1], maxn)
	if(n[4]!=maxn) nu  = rep(nu[1], maxn)
	if(n[5]!=maxn) tau = rep(tau[1], maxn)
	ans = double(maxn)
	sol = try(.C("c_st3_pst3", q = as.double(q), mu = as.double(mu), 
					sigma = as.double(sigma), nu = as.double(nu), 
					tau = as.double(tau), ans = ans, n = as.integer(maxn), 
					lower_tail = as.integer(lower.tail), 
					logr = as.integer(log.p), PACKAGE="gamlss.dist"), silent = TRUE)
	if(inherits(sol, 'try-error')){
		return(sol)
	} else{
		return(sol$ans)
	}
}
#-----------------------------------------------------------------  
qST3C <- function(p, mu=0, sigma=1, nu=1, tau=10, lower.tail = TRUE, log.p = FALSE)
{ 
	if (any(sigma < 0))  stop(paste("sigma must be positive", "\n", "")) 
	if (any(nu <= 0))  stop(paste("nu must be positive", "\n", ""))
	if (any(tau < 0))  stop(paste("tau must be positive", "\n", ""))  
	if (any(p <= 0)|any(p >= 1))  stop(paste("p must be between 0 and 1", "\n", ""))       
	n = c(length(p), length(mu), length(sigma), length(nu), length(tau))
	maxn = max(n)
	if(n[1]!=maxn) p = rep(p[1], maxn)
	if(n[2]!=maxn) mu    = rep(mu[1], maxn)
	if(n[3]!=maxn) sigma = rep(sigma[1], maxn)
	if(n[4]!=maxn) nu  = rep(nu[1], maxn)
	if(n[5]!=maxn) tau = rep(tau[1], maxn)
	ans = double(maxn)
	sol = try(.C("c_st3_qst3", p = as.double(p), mu = as.double(mu), 
					sigma = as.double(sigma), nu = as.double(nu), 
					tau = as.double(tau), ans = ans, n = as.integer(maxn), 
					lower_tail = as.integer(lower.tail), 
					logr = as.integer(log.p), PACKAGE="gamlss.dist"), silent = TRUE)
	if(inherits(sol, 'try-error')){
		return(sol)
	} else{
		return(sol$ans)
	}
}
#-----------------------------------------------------------------
# no benefit in using C code for now.
rST3C <- function(n, mu=0, sigma=1, nu=1, tau=10)
{
	if (any(sigma <= 0))  stop(paste("sigma must be positive", "\n", "")) 
	if (any(nu <= 0))  stop(paste("nu must be positive", "\n", ""))
	if (any(tau <= 0))  stop(paste("tau must be positive", "\n", ""))  
	if (any(n <= 0))  stop(paste("n must be a positive integer", "\n", ""))    
	n <- ceiling(n)
	p <- runif(n)
	r <- qST3(p,mu=mu,sigma=sigma,nu=nu,tau=tau)
	return(r)
}
