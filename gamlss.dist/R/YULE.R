#---------------------------#
#---- YULE DISTRIBUTION ----#
#---------------------------#

dYULE<-function (x, mu = 2, log = FALSE)
{
    if (any(mu < 0))
        stop(paste("mu must be > 0)", "\n", ""))
    if (any(x < 0))
        stop(paste("x must be >=0", "\n", ""))
    lx <- max(length(x), length(mu))
    mu<- rep(mu, length = lx)
    lambda <- (mu+1)/mu
    logfx <- lbeta(lambda+1, x+1) - lbeta(lambda, 1)
    if (log==FALSE) logfx <- exp(logfx)
    logfx
}

#Cumulative density function
pYULE<-function (q, mu = 2, lower.tail = TRUE, log.p = FALSE)
{
    if (any(mu < 0))
        stop(paste("mu must be > 0", "\n", ""))
    if (any(q < 0))
        stop(paste("q must be >=0", "\n", ""))
    ly <- max(length(q), length(mu))
    q <- rep(q, length = ly)
    mu<- rep(mu, length = ly)
    #s1<- seq(0, max(q))
    #cdf<- cumsum(dYUL(s1, mu=mu))
    #s2<-match(q,s1, nomatch=0)
    #cdf<- cdf[s2]
    cdf <- 1-((gamma(2+(1/mu))*gamma(2+q))/
           gamma(3+(1/mu)+q))
    if (lower.tail == TRUE) 
        cdf <- cdf
    else cdf = 1 - cdf
    if (log.p == TRUE) cdf <- -(lgamma(2+(1/mu))+lgamma(2+q) -
                              gamma(3+(1/mu)+q))
    cdf
}

#Quantile function
qYULE<-function (p, mu = 2, lower.tail = TRUE, log.p = FALSE, max.value = 10000)
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
                cumpro <- pYULE(j, mu= mu[i])
                QQQ[i] <- j
                if (p[i] <= cumpro)
                  break
            }
        }
    }
QQQ
}

#Random Generating Function
rYULE<- function(n, mu=2)
{
    if (any(mu < 0))
        stop(paste("mu must be > 0)", "\n", ""))
    if (any(n <= 0))
        stop(paste("n must be a positive integer", "\n", ""))
     n <- ceiling(n)
     p <- runif(n)
     r <- qYULE(p, mu=mu)
     r
}

#Distribution function
YULE<-function (mu.link = "log")
{
    mstats <- checklink(which.link="mu.link", 
    which.dist="Yule", link=substitute(mu.link),
     link.List="log")
    #One of these for each parameter to specify which link to use.

structure(list(family = c("YULE", "Yule"),      
        parameters = list(mu = TRUE),              
        nopar = 1,                        
        type = "Discrete",               
        mu.link = as.character(substitute(mu.link)),
        mu.linkfun = mstats$linkfun,
        mu.linkinv = mstats$linkinv,
        mu.dr = mstats$mu.eta,
        dldm = function(y, mu){
         lambda <- (mu+1)/mu
         dldm <- (digamma(lambda+1) - digamma(lambda+y+2)+(1/lambda))*(-1/(mu^2))
         dldm
          #browser()
        },                                     
        d2ldm2 = function(y, mu){
         # d2ldm2 <- 1/(mu*(mu-1))
           lambda <- (mu+1)/mu
           dldm <- (digamma(lambda+1) - digamma(lambda+y+2)+(1/lambda))*(-1/(mu^2))
           d2ldm2 <- -dldm^2
          d2ldm2 
        },
        G.dev.incr = function(y, mu, ...) 
            -2 * dYULE(y, mu = mu, log = TRUE),                
        rqres = expression(rqres(pfun = "pYULE", type = "Discrete", ymin = 0, y = y, mu = mu)),
        mu.initial = expression(mu <- rep(mean(y), length(y))),            
        mu.valid = function(mu) all(mu > 0) ,
        y.valid = function(y) all(y >=0)),        
        class = c("gamlss.family", "family"))
}