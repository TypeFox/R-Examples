#----------------------------------------------------------------------------------------
NET <- function (mu.link ="identity", sigma.link="log") 
{
    mstats <- checklink("mu.link",    "NET", substitute(   mu.link), c("inverse", "log", "identity", "own"))
    dstats <- checklink("sigma.link", "NET", substitute(sigma.link), c("inverse", "log", "identity", "own"))
    structure(
          list(family = c("NET", "Normal Exponential t"),
           parameters = list(mu=TRUE,sigma=TRUE, nu=FALSE, tau=FALSE), 
                nopar = 4, 
                 type = "Continuous",  
              mu.link = as.character(substitute(mu.link)), 
           sigma.link = as.character(substitute(sigma.link)), 
           mu.linkfun = mstats$linkfun, 
        sigma.linkfun = dstats$linkfun, 
           mu.linkinv = mstats$linkinv, 
        sigma.linkinv = dstats$linkinv,
                mu.dr = mstats$mu.eta, 
             sigma.dr = dstats$mu.eta, 
                 dldm = function(y,mu,sigma,nu,tau){ 
                                  k1 <- nu
                                  k2 <- tau
                                  tc <- (y-mu)/sigma
                               dldtc <- (abs(tc)<= k1)*(-tc) + 
                                        ((abs(tc)>k1)&(abs(tc)<= k2))*(-k1*sign(tc)) +
                                        (abs(tc) > k2)*(-k1*k2/tc)
                                dldm <- -dldtc/sigma
                                dldm 
                                    },
               d2ldm2 = function(y,mu,sigma,nu,tau) {
                                  k1 <- nu
                                  k2 <- tau
                                  c1 <- (1-2*pnorm(-k1))*sqrt(2*pi)
                                  c2 <- (2/k1)*exp(-((k1^2)/2))
                                  c3 <- 2*exp(-k1*k2+((k1^2)/2))/((k1*k2-1)*k1)
                                 ct  <- 1/(c1+c2+c3)
                                # tc  <- (y-mu)/sigma
                              d2ldm2 <- -ct*sqrt(2*pi)*(1-2*pnorm(-k1))+
                                        ((2*ct*k1)/(k1*k2+1))*exp(-k1*k2+(k1^2)/2)  
                              d2ldm2 <- d2ldm2/sigma^2},
                 dldd = function(y,mu,sigma,nu,tau) {
                                 k1 <- nu
                                 k2 <- tau
                                 tc <- (y-mu)/sigma
                              dldtc <- (abs(tc)<= k1)*(-tc) + 
                                        ((abs(tc)>k1)&(abs(tc)<= k2))*(-k1*sign(tc)) +
                                        (abs(tc) > k2)*(-k1*k2/tc)
                               dldd <- -(1+tc*dldtc)/sigma
                               dldd
                                    },
               d2ldd2 = function(y,mu,sigma,nu,tau) {
                                   k1 <- nu
                                   k2 <- tau
                                   c1 <- (1-2*pnorm(-k1))*sqrt(2*pi)
                                   c2 <- (2/k1)*exp(-((k1^2)/2))
                                   c3 <- 2*exp(-k1*k2+((k1^2)/2))/((k1*k2-1)*k1)
                                   ct <- 1/(c1+c2+c3)
                                  #tc  <- (y-mu)/sigma
                               d2ldd2 <- 2*ct*k1*exp(-(k1^2)/2)-
                                           ct*sqrt(2*pi)*(1-2*pnorm(-k1))+
                                           (2*ct*k1*(k2^2)/(k1*k2-1))*exp(-k1*k2+(k1^2)/2)
                               d2ldd2 <- (d2ldd2-1)/(sigma^2)
                               },
              d2ldmdd = function(y)  rep(0,length(y)),     
          G.dev.incr  = function(y,mu,sigma,nu,tau,...)  -2*dNET(y,mu,sigma,nu,tau, log=TRUE),  
                rqres = expression(rqres(pfun="pNET", type="Continuous", y=y, mu=mu, sigma=sigma, nu=nu, tau=tau)),
            mu.initial = expression(mu <- (y-mean(y))/2 ),
         sigma.initial = expression(sigma <- rep((sd(y)+0.001),length(y))),
            nu.initial = expression(nu <- rep(1.5,length(y))), # 1.5  default
           tau.initial = expression(tau <- rep(2,length(y))),  # 2   default
              mu.valid = function(mu) TRUE,  
           sigma.valid = function(sigma)  all(sigma > 0),
              nu.valid = function(nu) all(nu > 0), 
             tau.valid = function(tau, nu) all(tau > nu),
               y.valid = function(y) TRUE
          ),
            class = c("gamlss.family","family"))
}
#----------------------------------------------------------------------------------------
dNET <- function(x, mu=0, sigma=1, nu=1.5, tau=2, log=FALSE)
 {
          if (any(sigma <= 0))  stop(paste("sigma must be positive", "\n", "")) 
          if (any(nu <= 0))  stop(paste("nu must be positive", "\n", ""))  
          if (any(tau <= 0))  stop(paste("tau must be positive", "\n", ""))  
          if (any(tau < nu))  stop(paste(" tau must greater or equal than  nu", "\n", ""))  
        k1 <- nu
        k2 <- tau
        c1 <- (1-2*pnorm(-k1))*sqrt(2*pi)
        c2 <- (2/k1)*exp(-((k1^2)/2))
        c3 <- 2*exp(-k1*k2+((k1^2)/2))/((k1*k2-1)*k1)
       ct  <- 1/(c1+c2+c3)
       tc  <- (x-mu)/sigma
        d1 <- (abs(tc) <= k1)*(-(tc^2)/2)
        d2 <- ((abs(tc) > k1) & (abs(tc) <= k2))*(-k1*abs(tc)+((k1^2)/2))
        d3 <- (abs(tc) > k2)*(-k1*k2*log(abs(tc)/k2)-k1*k2+((k1^2)/2))
    loglik <- log(ct)-log(sigma)+d1+d2+d3 
      if(log==FALSE) ft  <- exp(loglik) else ft <- loglik 
       ft
  }      
#----------------------------------------------------------------------------------------
pNET <- function(q, mu=5, sigma=0.1, nu=1, tau=2)
 {  
          if (any(sigma <= 0))  stop(paste("sigma must be positive", "\n", "")) 
          if (any(nu <= 0))  stop(paste("nu must be positive", "\n", ""))  
          if (any(tau <= 0))  stop(paste("tau must be positive", "\n", ""))  
          if (any(tau < nu))  stop(paste(" tau must greater or equal than  nu", "\n", ""))  
        k1 <- nu[1]
        k2 <- tau[1]
        c1 <- (1-2*pnorm(-k1))*sqrt(2*pi)
        c2 <- (2/k1)*exp(-((k1^2)/2))
        c3 <- 2*exp(-k1*k2+((k1^2)/2))/((k1*k2-1)*k1)
        ct <- 1/(c1+c2+c3)
        tc <- (q-mu)/sigma
                                    #Fk2 is cdf up to the point -k2
       Fk2 <- (ct*k2^(k1*k2)/(k1*k2-1)) * exp(-k1*k2+(k1^2/2))* 
               abs(tc)^(-k1*k2+1) 
                                    #cf1 is cdf value up to the point -k2
      cdf1 <- (ct*k2/(k1*k2-1)) * exp(-k1*k2+(k1^2/2)) 
                                    #Fk1 is cdf up to the point -k1
       Fk1 <- cdf1 + ct * exp(k1^2/2) *(exp(-k1*abs(tc))-exp(-k1*k2))/k1
                                    #cf2 is cdf value up to the point -k1
      cdf2 <- cdf1 + ct * exp(k1^2/2) *(exp(-k1^2)-exp(-k1*k2))/k1 
                                    #F0 is cdf up to the point 0 
        F0 <- cdf2 + ct * sqrt(2*pi) * (pnorm(tc) - pnorm(-k1))
                                    #calclulate the cdf0=F(-tc) to get the cdf of
                                    #the positive tc as 1-F(-tc)
      cdf0 <- cdf2 + ct * sqrt(2*pi) * (pnorm(-tc) - pnorm(-k1))
                                    #cdf.tc is the cdf function of all the points
       cdf <- (tc <= -k2) * Fk2 + #negative tc
                ((tc > -k2) & (tc <= -k1)) * Fk1 + #negative tc
                ((tc > -k1) & (tc <= 0)) * F0 + #negative tc
                (tc > 0 ) * (1-cdf0) #positive tc
       cdf                                       
 }
#----------------------------------------------------------------------------------------
