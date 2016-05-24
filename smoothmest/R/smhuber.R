smoothm <- function(y, method="smhuber",
                    k=0.862,
                    sn=sqrt(2.046/length(y)), tol = 1e-06, s=mad(y),
                    init="median")
{
  estim <- switch(method,
                  huber=sehuber(y,k,tol,s,init),
                  smhuber=smhuber(y,k,sn,tol,s,init=init),
                  bisquare=mbisquare(y,k,tol,s,init),
                  smbisquare=smbisquare(y,k,tol,sn,s,init),
                  cauchy=mlcauchy(y,tol,s),
                  smcauchy=smcauchy(y,tol,sn,s),
                  smmed=smhuber(y,k,sn,tol,s,smmed=TRUE,init))
  out <- list(mu=estim$mu,method=method,k=k,sn=sn,tol=tol,s=s)
  out
}

sehuber <- function (y, k = 0.862, tol = 1e-06, s=mad(y), init="median") 
{
    y <- y[!is.na(y)]
    n <- length(y)
    if (init=="median")
      mu <- median(y)
    if (init=="mean")
      mu <- mean(y)
    repeat {
        yy <- pmin(pmax(mu - k * s, y), mu + k * s)
        mu1 <- sum(yy)/n
        if (abs(mu - mu1) < tol * s) 
            break
        mu <- mu1
    }
    list(mu = mu, s = s)
}

psicauchy <- function(x) 2*x/(1+x^2)
psidcauchy <- function(x) (2*(1+x^2)-4*x^2)/(1+x^2)^2
# avcauchy <- function(x) psicauchy(x)^2*dcauchy(x)
likcauchy <- function(x,mu) prod(dcauchy(x-mu))
flikcauchy <- function(y,x,mu,sn){
  out <- c()
  for (i in 1:length(y))
    out[i] <- sum(log(dcauchy(x-y[i]-mu)))*dnorm(y[i],sd=sn)
  out
}

smtfcauchy <- function(x,mu,sn)
  integrate(flikcauchy, -Inf, Inf, x=x, mu=mu, sn=sn)$value

smcipsi <- function(y, x, sn=sqrt(2/length(x)))
  psicauchy(x-y)*dnorm(y,sd=sn)

smcipsid <- function(y, x, sn=sqrt(2/length(x)))
  psidcauchy(x-y)*dnorm(y,sd=sn)

smcpsi <- function(x, sn=sqrt(2/length(x))){
  out <- c()
  for (y in x)
    out <- c(out,integrate(smcipsi, -Inf, Inf, x=y, sn=sn)$value)
  out
}

smcpsid <- function(x, sn=sqrt(2/length(x))){
  out <- c()
  for (y in x)
    out <- c(out,integrate(smcipsid, -Inf, Inf, x=y, sn=sn)$value)
  out
}




smbpsi <- function(y, x, k=4.685, sn=sqrt(2/length(x))){
  require(MASS)
  (x-y)*psi.bisquare(x-y,c=k)*dnorm(y,sd=sn)
}

smbpsid <- function(y, x, k=4.685, sn=sqrt(2/length(x))){
  require(MASS)
  psi.bisquare(x-y,c=k,deriv=1)*dnorm(y,sd=sn)
}

smbpsii <- function(x, k=4.685, sn=sqrt(2/length(x))){
  out <- c()
  for (y in x)
    out <- c(out,integrate(smbpsi, -Inf, Inf, x=y, k=k, sn=sn)$value)
  out
}

smbpsidi <- function(x, k=4.685, sn=sqrt(2/length(x))){
  out <- c()
  for (y in x)
    out <- c(out,integrate(smbpsid, -Inf, Inf, x=y, k=k, sn=sn)$value)
  out
}

mbisquare <- function (y, k=4.685, tol = 1e-06, s=mad(y), init="median") 
{
    require(MASS)
    y <- y[!is.na(y)]
    n <- length(y)
    if (init=="median")
      mu <- median(y)
    if (init=="mean")
      mu <- mean(y)
    ic <- 0
    repeat {
        ic <- ic+1
        s1 <- sum((y-mu)*psi.bisquare((y-mu)/s,c=k))
        s2 <- sum(psi.bisquare((y-mu)/s,c=k,deriv=1))
        if (abs(s2)>tol)
          mu1 <- mu+s1/s2
        else
          mu1 <- mu
        if (mu1>max(y) | mu1<min(y) | ic>200)
          break
        if (abs(mu - mu1) < tol * s) 
            break
        mu <- mu1
    }
    list(mu = mu, s = s)
}

smbisquare <- function (y, k=4.685, tol = 1e-06, sn=sqrt(1.0526/length(y)), s=mad(y), init="median") 
{
    y <- y[!is.na(y)]
    n <- length(y)
    if (init=="median")
      mu <- median(y)
    if (init=="mean")
      mu <- mean(y)
    ic <- 0
    repeat {
        ic <- ic+1
        s1 <- sum(s*smbpsii((y-mu)/s, k=k, sn=sn))
        s2 <- sum(smbpsidi((y-mu)/s, k=k, sn=sn))
        if (abs(s2)>tol)
          mu1 <- mu+s1/s2
        else
          mu1 <- mu
        if (mu1>max(y) | mu1<min(y) | ic>200)
          break
        if (abs(mu - mu1) < tol * s) 
            break
        mu <- mu1
    }
    list(mu = mu, s = s)
}



mlcauchy <- function (y, tol = 1e-06, s=mad(y)) 
{
    y <- y[!is.na(y)]
    n <- length(y)
    me <- mu <- median(y)
    ic <- 0
    lmed <- likcauchy(y,mu)
    repeat {
        ic <- ic+1
        s1 <- sum(s*psicauchy((y-mu)/s))
        s2 <- sum(psidcauchy((y-mu)/s))
        mu1 <- mu+s1/s2
        if (mu1>max(y) | mu1<min(y) | ic>200)
          break
        if (abs(mu - mu1) < tol * s) 
            break
        mu <- mu1
    }
    lmu <- likcauchy(y,mu)
    if (lmu<=lmed){
      md <- abs(rank(y)-(n+1)/2)
      bi <- 0
      md[md==0] <- n
      mus <- me
      for (q in (1+2*(0:(ceiling(n/2)-2)))){
        for (j in 0:1){
          i <- order(md)[q+j]
          mu <- y[i]
          lmunew <- likcauchy(y,mu)
          if (lmunew>lmed){
            bi <- i
            lmed <- lmunew
          }
          ic <- 0
          repeat {
            ic <- ic+1
            s1 <- sum(s*psicauchy((y-mu)/s))
            s2 <- sum(psidcauchy((y-mu)/s))
            mu1 <- mu+s1/s2
            if (mu1>max(y) | mu1<min(y) | ic>200)
              break
            if (abs(mu - mu1) < tol * s) 
              break
            mu <- mu1
          }
          lk <- likcauchy(y,mu)
          if(lk>lmu){
            lmu <- lk
            mus <- mu
          }
        }
        mu <- mus
        if (lmu>lmed)
          break
      }
      if (lmu<=lmed){
        mu <- max(y)
        mu1 <- min(y)
        sa <- sum(psicauchy((y-mu)/s))
        sb <- sum(psicauchy((y-mu1)/s))
        repeat{
          mu2 <- (mu+mu1)/2
          sc <- sum(psicauchy((y-mu2)/s))
          if (sign(sc)==sign(sa)) mu <- mu2
          else mu1 <- mu2
          if (abs(mu - mu1) < tol * s) 
            break
        }
        lmu <- likcauchy(y,mu)
      }
    }
    if (lmu>lmed) mu1 <- mu
    else mu1 <- ifelse(bi==0,me,y[bi])
    list(mu = mu1, s = s)
}

smcauchy <- function (y, tol = 1e-06, sn=sqrt(2/length(y)), s=mad(y)) 
{
    y <- y[!is.na(y)]
    n <- length(y)
    me <- mu <- median(y)
    ic <- 0
    lmed <- smtfcauchy(y,mu,sn)
    repeat {
        ic <- ic+1
        mu1 <- mu+sum(s*smcpsi((y-mu)/s, sn=sn))/
                  sum(smcpsid((y-mu)/s, sn=sn))
        if (mu1>max(y) | mu1<min(y) | ic>200)
          break
        if (abs(mu - mu1) < tol * s) 
            break
        mu <- mu1
    }
    lmu <- smtfcauchy(y,mu,sn)
    if (lmu<=lmed){
      md <- abs(rank(y)-(n+1)/2)
      bi <- 0
      md[md==0] <- n
      mus <- me
      for (q in (1+2*(0:(ceiling(n/2)-2)))){
        for (j in 0:1){
          i <- order(md)[q+j]
          mu <- y[i]
          lmunew <- smtfcauchy(y,mu,sn)
          if (lmunew>lmed){
            bi <- i
            lmed <- lmunew
          }
          ic <- 0
          repeat {
            ic <- ic+1
            mu1 <- mu+sum(s*smcpsi((y-mu)/s, sn=sn))/
                  sum(smcpsid((y-mu)/s, sn=sn))
            if (mu1>max(y) | mu1<min(y) | ic>200)
              break
            if (abs(mu - mu1) < tol * s) 
              break
            mu <- mu1
          }
          lk <- smtfcauchy(y,mu,sn)
          if(lk>lmu){
            lmu <- lk
            mus <- mu
          }
        }
        mu <- mus
        if (lmu>lmed)
          break
      }
    }
    if (lmu>lmed) mu1 <- mu
    else mu1 <- ifelse(bi==0,me,y[bi])
    list(mu = mu1, s = s)
}



smpsi <- function(x,k=0.862,sn=sqrt(2/length(x)))
  {
    k*pnorm(x,k,sn)-k*(1-pnorm(x,-k,sn))+x*(pnorm(x,-k,sn)-pnorm(x,k,sn))+
                                            sn*(dnorm((x+k)/sn)-dnorm((x-k)/sn))
  }

smpmed <- function(x,sn=sqrt(1/5))
{
  -pnorm(0,x,sn)+pnorm(0,-x,sn)
}

# For median=TRUE:
# Double exp. distribution with 1.4826MAD=1, lambda=1/(1.4826*log(2)),
# V=1/(4*(1/2lambda)^2)=1.056

smhuber <- function (y, k = 0.862, sn=sqrt(2.046/length(y)), tol = 1e-06, s=mad(y),
                     smmed=FALSE, init="median") 
{
    y <- y[!is.na(y)]
    n <- length(y)
    if (init=="median")
      mu <- median(y)
    if (init=="mean")
      mu <- mean(y)
    repeat {
        if(smmed){
          w <- smpmed(y-mu,sn=sn)/(y-mu)
          w[abs(w)==Inf] <- NA
          w[is.na(w)] <- max(w,na.rm=TRUE)
          w[is.na(w)] <- 1
        }          
        else{  
          w <- smpsi(y-mu,k=k,sn=sn)/(y-mu)
          w[abs(w)==Inf] <- NA
          w[is.na(w)] <- max(w,na.rm=TRUE)
          w[is.na(w)] <- 1
        }          
        yy <- w*y
        mu1 <- sum(yy)/sum(w)
        if (abs(mu - mu1) < tol * s) 
            break
        mu <- mu1
    }
    list(mu = mu, s = s)
}

dens <- function(x, dfunction, ...){
  prod(dfunction(x, ...))
}
pdens <- function(z, x, dfunction, ...)
{
  out <- c()
  for (j in 1:length(z)){
    out[j] <- z[j]*dens(x-z[j], dfunction, ...)
  }
  out
}

sdens <- function(z, x, dfunction, ...) 
{
  out <- c()
  for (j in 1:length(z))
    out[j] <- dens(x-z[j], dfunction, ...)
  out
}

pitman <- function(y, d=ddoublex, lower=-Inf, upper=Inf, s=mad(y), ...)
{
  n <- length(y)
  z <- y/s
  out <- s*integrate(pdens, lower, upper, x=z, dfunction=d, ...)$value/
         integrate(sdens, lower, upper, x=z, dfunction=d, ...)$value
  if (is.na(out)) out <- median(y)
  out
}
  
edhuber <- function(x, k=0.862, mu=0, sigma=1)
{
  z <- (x-mu)/sigma
  fk <- dnorm(k)
  eps <- 1-1/(pnorm(k)-pnorm(-k)+2*fk/k)
  val=c()
  for (y in z){
    if (y<(-k))
      val <- c(val,(1-eps)*fk*exp(k*(y+k)))
    else{
      if (y>k)
        val <- c(val,(1-eps)*fk*exp(-k*(y-k)))
      else
        val <- c(val,(1-eps)*dnorm(y))
    }
  }
  out <- list(val=val/sigma,eps=eps)
  out
}

dhuber <- function(x, k=0.862, mu=0, sigma=1) edhuber(x, k, mu, sigma)$val

rhuber <- function(n,k=0.862, mu=0, sigma=1)
{
  l <- c()
  while(length(l)<n){
    x <- rexp(1)
    s <- sample(c(-1,1),size=1)
    y <- s*x/k
    u <- runif(1)
    if (abs(y)>=k | u<=exp(k*abs(y)-(k*k+y*y)/2))
      l <- c(l,y)
  }
  sigma*l+mu
}

ddoublex <- function(x, mu=0, lambda=1)
  exp(-abs(x-mu)/lambda)/(2*lambda)

rdoublex <- function(n,mu=0,lambda=1)
{
    x <- rexp(n)
    s <- sample(c(-1,1),size=n,replace=TRUE)
    y <- s*x
    lambda*y+mu
}
  


