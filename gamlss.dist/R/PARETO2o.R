######################################
##### PARETO TYPE 2 DISTRIBUTION #####
######################################
#-------------------------------------------------------------------------------
# Robert Rigby, Fiona McElduf, Mikis Stasinopoulos, Vlasios Voudouris 
#################################################################################
#-------------------------------------------------------------------------------
#Probability distribution function
dPARETO2o <- function(x, mu = 1, sigma = 0.5, log = FALSE)
{
   if (any(mu < 0)) stop(paste("mu must be positive", "\n", "")) 
   if (any(sigma <= 0))   stop(paste("sigma must be positive", "\n", ""))  
   if (any(x < 0)) stop(paste("x must be greater than 0", "\n", ""))
   lfy <- log(sigma) + sigma*log(mu) - (sigma+1)*log(x+mu)
   if (log == FALSE) fy <- exp(lfy) else fy <- lfy
   fy
}
#--------------------------------------------------------------------------------
#Cumulative density function
pPARETO2o <- function(q, mu = 1, sigma = 0.5, lower.tail = TRUE, log.p = FALSE)
{
   if (any(mu <= 0)) stop(paste("mu must be positive", "\n", "")) 
   if (any(sigma <= 0)) stop(paste("tau must be positive", "\n", ""))           
   if (any(q < 0)) stop(paste("q must be be greater than 0", "\n", ""))   
   cdf <- 1 - ((mu/(mu+q))^(sigma))
   if (lower.tail == TRUE) cdf <- cdf  
   else cdf <- 1 - cdf
   if (log.p == FALSE) cdf <- cdf
   else cdf < - log(cdf)
   cdf
}   
#-------------------------------
#Quantile-inverse cdf  
qPARETO2o <- function(p, mu = 1, sigma = 0.5, lower.tail = TRUE, log.p = FALSE)
{
  if (any(mu < 0))  stop(paste("mu must be positive", "\n", "")) 
#    if (any(nu < 0))  stop(paste("nu must be positive", "\n", ""))  
    if (any(sigma < 0))  stop(paste("sigma must be positive", "\n", ""))  
    if (log.p==TRUE) p <- exp(p) else p <- p
    if (any(p <= 0)|any(p >= 1))  stop(paste("p must be between 0 and 1", "\n", ""))       
    if (lower.tail==TRUE) p <- p else p <- 1-p
   # w <- qf(p,2,2/sigma)   
   #  q1 <- mu*(((sigma)*w))
    q <-  mu*((1-p)^(-(1/sigma))-1)   
    q
}
#--------------------------------------------------------------------------------
#Random generation 
rPARETO2o <- function(n, mu = 1, sigma = 0.5)
{
   if (any(mu <= 0)) stop(paste("mu must be positive", "\n", "")) 
   if (any(sigma <= 0)) stop(paste("sigma must be positive", "\n", "")) 
   if (any(n <= 0)) stop(paste("n must be a positive integer", "\n", ""))  
   n <- ceiling(n)
   p <- runif(n)
   r <- qPARETO2o(p, mu = mu, sigma = sigma)
   r 
}
#-------------------------------------------------------------------------------
#Gamlss Family Function
PARETO2o <- function (mu.link = "log", sigma.link = "log") 
{
    mstats <- checklink("mu.link", "Pareto Type 2", substitute(mu.link), 
        c("inverse", "log", "identity", "own"))
    dstats <- checklink("sigma.link", "Pareto Type 2", substitute(sigma.link), 
        c("inverse", "log", "identity", "own"))
    structure(
       list(family = c("PARETO2o", "Pareto Type 2"), 
        parameters = list(mu = TRUE, sigma = TRUE), 
             nopar = 2, 
              type = "Continuous", 
           mu.link = as.character(substitute(mu.link)), 
        sigma.link = as.character(substitute(sigma.link)), 
        mu.linkfun = mstats$linkfun, 
     sigma.linkfun = dstats$linkfun, 
        mu.linkinv = mstats$linkinv, 
     sigma.linkinv = dstats$linkinv, 
             mu.dr = mstats$mu.eta, 
          sigma.dr = dstats$mu.eta, 
              dldm = function(y, mu, sigma)
              {
               dldm <- (sigma/mu)-((sigma+1)/(y+mu))
               dldm
                },      
            d2ldm2 = function(y, mu, sigma)
              {
              #d2ldm2 <- -(sigma/(mu^2*(sigma+2)))
              dldm <- (sigma/mu)-((sigma+1)/(y+mu))
              d2ldm2 <- -dldm^2
              d2ldm2
              },
            dldd = function(y, mu, sigma)  
              {
              dldd <- (1/sigma)+log(mu)-log(y+mu)
              dldd
              },
            d2ldd2 = function(y, mu, sigma)      
              {
              #d2ldd2 <- -(1/sigma^2) 
              dldd <- (1/sigma)+log(mu)-log(y+mu)
              d2ldd2 <- -dldd^2
              d2ldd2  
              },
            d2ldmdd = function(y, mu, sigma) #
              {
              #d2ldmdd <- -1/(mu+mu*sigma)
              dldm <- (sigma/mu)-((sigma+1)/(y+mu))
              dldd <- (1/sigma)+log(mu)-log(y+mu)
              d2ldmdd <- -dldm*dldd
              d2ldmdd 
              },
        G.dev.incr = function(y, mu, sigma, ...) -2 * 
            dPARETO2o(y, mu, sigma, log = TRUE), 
        rqres = expression(rqres(pfun = "pPARETO2o", 
            type = "Continuous", y = y, mu = mu, sigma = sigma)), 
        mu.initial = expression({mu <- rep(mean(y), length(y))}), 
        sigma.initial = expression({sigma <- rep(sd(y), length(y))}), 
        mu.valid = function(mu) all(mu > 0), 
        sigma.valid = function(sigma) all(sigma > 0), 
        y.valid = function(y) TRUE), 
        class = c("gamlss.family", "family"))
}
