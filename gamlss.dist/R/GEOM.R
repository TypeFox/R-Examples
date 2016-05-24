#--------------------------------#
#---- GEOMETRIC DISTRIBUTION ----#
#--------------------------------#

dGEOM<-function (x, mu = 2, log = FALSE)
{
    if (any(mu < 0))
        stop(paste("mu must be > 0)", "\n", ""))
    if (any(x < 0))
        stop(paste("x must be >=0", "\n", ""))
   # lx <- max(length(x), length(mu))
   # mu<- rep(mu, length = lx)
    #prob <- 
    #if (length(x)<100) browser()
    #logfx <- x*log(1-(1/(mu+1))) + log(1/(mu+1))
    logfx <- x*log(mu/(mu+1)) + log(1/(mu+1))
   #  cat("x", x, "\n")
   #  browser()
    if (log==FALSE) logfx <- exp(logfx)
    logfx
}

#Cumulative density function
pGEOM<-function (q, mu = 2, lower.tail = TRUE, log.p = FALSE)
{
    if (any(mu < 0))
        stop(paste("mu must be > 0", "\n", ""))
    if (any(q < 0))
        stop(paste("q must be >=0", "\n", ""))
   # ly <- max(length(q), length(mu))
  #  q <- rep(q, length = ly)
  #  mu<- rep(mu, length = ly)
  #  s1<- seq(0, max(q))
  #  cdf<- cumsum(dGEOM(s1, mu=mu))
  #  s2<-match(q,s1, nomatch=0)
   # cdf<- cdf[s2]
     cdf <- 1-(mu/(mu+1))^(q+1)  
     cdf <- if (log.p==FALSE) cdf else log(cdf)
     cdf <- if (lower.tail == TRUE) cdf   else  1 - cdf
    if (log.p == TRUE) cdf<-log(cdf)
    cdf
}

#Quantile function
qGEOM<-function (p, mu = 2, lower.tail = TRUE, log.p = FALSE, max.value = 10000)
{
    if (any(p < 0) | any(p > 1.0001))
        stop(paste("p must be in [0,1]", "\n", ""))
    if (any(mu < 0))
        stop(paste("mu must be > 0)", "\n", ""))     
if (lower.tail) p <- p
else p <- 1 - p
    ly <- max(length(p), length(mu))
    p <- rep(p, length = ly)
    QQQ <- rep(0, length = ly)
    mu <- rep(mu, length = ly)
    for (i in seq(along = p)) {
        cumpro <- 0
        if (p[i] + 1e-09 >= 1)
            QQQ[i] <- Inf
        else {
            for (j in seq(from = 0, to = max.value)) {
                cumpro <- pGEOM(j, mu= mu[i])
                QQQ[i] <- j
                if (p[i] <= cumpro)
                  break
            }
        }
    }
QQQ
}

#Random Generating Function
rGEOM<- function(n, mu=2)
{
    if (any(mu < 0))
        stop(paste("mu must be > 0)", "\n", ""))
    if (any(n <= 0))
        stop(paste("n must be a positive integer", "\n", ""))
     n <- ceiling(n)
     p <- runif(n)
     r <- qGEOM(p, mu=mu)
     r
}

#Distribution function
GEOM<-function (mu.link = "log")
{
mstats <- checklink("mu.link", "Geometric", substitute(mu.link), 
        c("log", "probit", "cloglog", "cauchit", "log", "own"))
structure(list(family = c("GEOM", "Geometric"),      
        parameters = list(mu = TRUE),     
        nopar = 1,                       
        type = "Discrete",               
        mu.link = as.character(substitute(mu.link)),  
        mu.linkfun = mstats$linkfun,
        mu.linkinv = mstats$linkinv,
        mu.dr = mstats$mu.eta,
        dldm = function(y, mu){
          dldm <- (y - mu)/(mu + (mu^2))
          dldm
        },                                             
        d2ldm2 = function(mu){
          d2ldm2 <- -1/(mu+(mu^2))
          d2ldm2
        }, 
        G.dev.incr = function(y, mu, ...) -2*dGEOM(x=y, mu=mu, log = TRUE),                 
        rqres = expression(rqres(pfun = "pGEOM", type = "Discrete", 
                            ymin = 0, y = y, mu = mu)),   
        mu.initial = expression(mu <- rep(mean(y), length(y))),            
        mu.valid = function(mu) all(mu > 0) , 
        y.valid = function(y) all(y >=0)),      
        class = c("gamlss.family", "family"))
}

