# Density functions of SSMN Distributions
#require(moments)
#require(truncdist)

dssmn <- function(x, location=0, scale=1, shape=0, nu= 1, gama=1, dp=NULL, family="sn")
{
 if(family=="sn")
 {
  if(!is.null(dp))
  {
   if(!missing(shape)){stop("You cannot set both component parameters and dp")}
   location <- dp[1]
   scale    <- dp[2]
   shape    <- dp[3]
  }
  if(scale<=0){stop("Parameter scale must be positive")}
  if(any(is.na(x))){x <- x[-which(is.na(x))]}
  z <- (x-location)/sqrt(scale)
  y <- 2*dnorm(z)*pnorm(shape*z)/sqrt(scale)
 }

 if(family=="stn")
 {
  if(!is.null(dp))
  {
   if(!missing(shape)){stop("You cannot set both component parameters and dp")}
   location <- dp[1]
   scale    <- dp[2]
   shape    <- dp[3]
   nu       <- dp[4]
  }
  if (nu<=0){stop("Parameter nu must be positive")}
  if(scale<=0){stop("Parameter scale must be positive")}
  if (any(is.na(x))){ x <- x[-which(is.na(x))] }
  z    <- (x-location)/sqrt(scale)
  cte  <- gamma((nu+1)/2)/gamma(nu/2)/sqrt(nu*pi)
  pdft <- cte*(1+z^2/nu)^(-(nu+1)/2)
  y    <- 2*pdft*pnorm(shape*z)/sqrt(scale)
 }

 if(family=="stn")
 {
  if(!is.null(dp))
  {
   if(!missing(shape)){stop("You cannot set both component parameters and dp")}
   location <- dp[1]
   scale    <- dp[2]
   shape    <- dp[3]
   nu       <- dp[4]
  }
  if(nu<=0){stop("Parameter nu must be positive")}
  if(scale<=0){stop("Parameter scale must be positive")}
  if(any(is.na(x))){x  <- x[-which(is.na(x))]}
  z    <- (x-location)/sqrt(scale)
  n    <- length(x)
  cte  <- nu/((2*pi*scale)^(1/2))
  d    <- z^2
  fdps <- cte*gamma(nu+1/2)*pgamma(1,nu+1/2,scale=2/d)/((d/2)^(nu+1/2))
  fdps[which(z==0)] <- cte/(nu+1/2)
  y  <- 2*fdps*pnorm(shape*z)
 }

 if(family=="stn")
 {
  if(!is.null(dp))
  {
   if(!missing(shape)){stop("You cannot set both component parameters and dp")}
   location <- dp[1]
   scale    <- dp[2]
   shape    <- dp[3]
   nu       <- dp[4]
   gama     <- dp[5]
  }
  if(nu<=0|nu>=1){stop("Parameter nu must be between 0 and 1.0")}
  if(gama<=0|gama>=1){stop("Parameter gama must be between 0 and 1.0")}
  if(scale<=0){stop("Parameter scale must be positive")}
  if(any(is.na(x))){x <- x[-which(is.na(x))]}
  z  <- (x-location)/sqrt(scale)
  z2 <- z*sqrt(gama)
  y  <- 2*(nu*dnorm(x, mean=location, sd=sqrt(scale/gama))+(1-nu)*dnorm(x,mean=location,sd=sqrt(scale)))*pnorm(shape*z)
 }

 if(family=="sep")
 {
  if(!is.null(dp))
  {
   if(!missing(shape)){stop("You cannot set both component parameters and dp")}
   location <- dp[1]
   scale    <- dp[2]
   shape    <- dp[3]
   nu       <- dp[4]
  }
  if(nu<=0.5|nu>1){stop("Parameter nu must be 0.5 < nu <= 1.0")}
  if(scale<=0){stop("Parameter scale must be positive")}
  if(any(is.na(x))){x<- x[-which(is.na(x))]}
  nu1   <- 1/(2*nu)
  z     <- (x-location)/sqrt(scale)
  d     <- z^2
  cte   <- nu/(2^nu1)/gamma(nu1)/sqrt(scale)
  pdfep <- cte*exp(-0.5*d^(nu))
  y     <- 2*pdfep*pnorm(shape*z)
 }

 return(list(y=y,family=family))
}


###############################################################################
# Cumulative Distributions functions of SSMN Distributions

pssmn <- function(q, location=0, scale=1, shape=0, nu= 1, gama=1, dp=NULL, family="sn")
{
 x    <- q
 nrep <- length(x)
 if(family=="sn")
 {
  if(!is.null(dp))
  {
   if(!missing(shape)){stop("You cannot set both component parameters and dp")}
   location <- dp[1]
   scale    <- dp[2]
   shape    <- dp[3]
  }else{
   if(scale<=0){stop("Parameter scale must be positive")}
   delta <- shape/sqrt(1+shape^2)
   z     <- (x-location)/sqrt(scale)
   p     <- numeric(length(x))
   for (i in 1:length(x))
   {
    if(abs(x[i]) == Inf){ p[i] <- (1+sign(x[i]))/2
    }else{
     covi <- matrix(c(1,-delta,-delta,1),2,2)
     Phi2 <- pmnorm(c(z[i],0),c(0,0),covi)
     p[i] <- 2*Phi2
         }
   }
   y <- pmax(0,pmin(1,p))
      }
 }

 if(family=="stn")
 {
  if(!is.null(dp))
  {
   if(!missing(shape)){stop("You cannot set both component parameters and dp")}
    location <- dp[1]
    scale    <- dp[2]
    shape    <- dp[3]
    nu       <- dp[4]
  }else{
   if(nu<=0){stop("Parameter nu must be positive")}
   if(scale<=0){stop("Parameter scale must be positive")}
   delta <- shape/sqrt(1+shape^2)
   U     <- rgamma(nrep,shape=nu/2,scale=2/nu)
   z     <- (x-location)/sqrt(scale)
   p     <- numeric(length(x))
   for(i in 1:length(x))
   {
    if(abs(x[i]) == Inf){ p[i] <- (1+sign(x[i]))/2
    }else{
      Phi2  <- numeric(nrep)
      for (j in 1:nrep)
      {
       kappaj  <- 1/U[j]
       fj      <- (1+shape^2*kappaj)/((1+shape^2)*kappaj)
       covi    <- kappaj*matrix(c(1,-delta,-delta,fj),2,2)
       Phi2[j] <- pmnorm(c(z[i],0),c(0,0),covi)
       }
       p[i]      <- 2*mean(Phi2)
        }
   }
    y <- pmax(0,pmin(1,p))
     }
 }

 if(family=="ssl")
 {
  if(!is.null(dp))
  {
   if(!missing(shape)){stop("You cannot set both component parameters and dp")}
   location <- dp[1]
   scale    <- dp[2]
   shape    <- dp[3]
   nu       <- dp[4]
  }else{
    if(nu<=0){stop("Parameter nu must be positive")}
    if(scale<=0){stop("Parameter scale must be positive")}
    delta <- shape/sqrt(1+shape^2)
    V     <- runif(nrep)
    U     <- V^(1/nu)
    z     <- (x-location)/sqrt(scale)
    p     <- numeric(length(x))
    for (i in 1:length(x))
    {
     if(abs(x[i]) == Inf) {p[i] <- (1+sign(x[i]))/2
     }else{
       Phi2 <- numeric(nrep)
       for (j in 1:nrep)
       {
        kappaj  <- 1/U[j]
        fj      <- (1+shape^2*kappaj)/((1+shape^2)*kappaj)
        covi    <- kappaj*matrix(c(1,-delta,-delta,fj),2,2)
        Phi2[j] <- pmnorm(c(z[i],0),c(0,0),covi)
       }
        p[i]    <- 2*mean(Phi2)
         }
    }
      y <- pmax(0,pmin(1,p))
      }
 }

 if(family=="scn")
 {
  if(!is.null(dp))
  {
   if(!missing(shape)) {stop("You cannot set both component parameters and dp")}
   location <- dp[1]
   scale    <- dp[2]
   shape    <- dp[3]
   nu       <- dp[4]
  }else{
    if(nu<=0){stop("Parameter nu must be positive")}
    if (scale<=0){stop("Parameter scale must be positive")}
    delta <- shape/sqrt(1+shape^2)
    V     <- rbinom(nrep,1,nu)
    U     <- gama*V+1-V
    z     <- (x-location)/sqrt(scale)
    p     <- numeric(length(x))
    for (i in 1:length(x))
    {
     if(abs(x[i]) == Inf) {p[i] <- (1+sign(x[i]))/2
     }else{
       Phi2 <- numeric(nrep)
       for (j in 1:nrep)
       {
         kappaj  <- 1/U[j]
         fj      <- (1+shape^2*kappaj)/((1+shape^2)*kappaj)
         covi    <- kappaj*matrix(c(1,-delta,-delta,fj),2,2)
         Phi2[j] <- pmnorm(c(z[i],0),c(0,0),covi)
       }
         p[i]    <- 2*mean(Phi2)
          }
    }
      y <- pmax(0,pmin(1,p))
      }
 }

 if(family=="sep") #tem um problema no codigo p[i]  <-
 {
  if(!is.null(dp))
  {
   if(!missing(shape)){stop("You cannot set both component parameters and dp")}
   location <- dp[1]
   scale    <- dp[2]
   shape    <- dp[3]
   nu       <- dp[4]
  }else{
    if(nu<=0.5|nu>1){stop("Parameter nu must be 0.5 < nu <= 1.0")}
    if(scale<=0){stop("Parameter scale must be positive")
      }
      #   z <- (x-location)/sqrt(scale)
      p <- numeric(length(x))
      for (i in 1:length(x)){
        p[i]  <-
          if(abs(x[i]) == Inf)   (1+sign(x[i]))/2
        else{
          integrate(dsep, -Inf, x[i], location=location, scale=scale, shape = shape, nu = nu)$value
        }
      }
      pmax(0,pmin(1,p))
    }
  }
}

dsep <- function(x, location=0, scale=1, shape=0, nu=1, dp=NULL)
{
  if(!is.null(dp)) {
    if(!missing(shape))
      stop("You cannot set both component parameters and dp")
    location <- dp[1]
    scale <- dp[2]
    shape <- dp[3]
    nu<- dp[4]
  }
  if (nu<0.5|nu>1) {
    stop("Parameter nu must be between 0.5 and 1.0")
  }
  if(scale<0) {
    stop("Parameter scale must be positive")
  }
  if (any(is.na(x))) x<- x[-which(is.na(x))]
  nu1<-1/(2*nu)
  z<-(x-location)/sqrt(scale)
  d<-z^2
  cte<-nu/(2^nu1)/gamma(nu1)/sqrt(scale)
  yt<-cte*exp(-0.5*d^(nu))
  y<-2*yt*pnorm(shape*z)
  return(y)
}

#pssmn(1, location=0, scale=1, shape=0, nu= 1, gama=1, dp=NULL, family="sn")
##############################################################################
# Quantile functions of SSMN Distributions #Revisar

qssmn <- function(p, location=0, scale=1, shape=0, nu= 1, gama=1, dp=NULL, family="sn")
{
 if(family=="sn")
 {
  if(!is.null(dp))
  {
   if(!missing(shape)){stop("You cannot set both component parameters and dp")}
   location <- dp[1]
   scale    <- dp[2]
   shape    <- dp[3]
  }else{
    if(scale<=0){stop("Parameter scale must be positive")}
    na   <- is.na(p) | (p < 0) | (p > 1)
    zero <- (p == 0)
    one  <- (p == 1)
    q    <- numeric(length(p))
    aux1 <- which(na | one | zero)
    if(length(p)!=length(aux1))
    {
     aux2 <- c(1:length(p))
     # aux2<- ifelse(length(aux1)!=0,aux2[-c(aux1)],aux2)
     for (i in aux2)
     {
          q[i]  <-
            if (shape==0) qnorm(p[i])
          else         optimize(function(x) abs(psnn(x,0,1,shape)-p[i]),c((1.125 - 2.575*sign(shape))*sign(shape),(1.125 + 2.575*sign(shape))*sign(shape)))$minimum
     }
    }
        }
    q <- replace(q, na, NA)
    q <- replace(q, zero, -Inf)
    q <- replace(q, one, Inf)
    return(location+sqrt(scale)*q)
  }


qstn <- function (p,location=0, scale=1, shape=0, nu=30, dp=NULL) {
    if(!is.null(dp)) {
      if(!missing(shape))
        stop("You cannot set both component parameters and dp")
      location <- dp[1]
      scale <- dp[2]
      shape <- dp[3]
      nu <- dp[4]
    }
    else
    {
      if (nu<=0) {
        stop("Parameter nu must be positive")
      }
      if (scale<=0) {
        stop("Parameter scale must be positive")
      }
      na   <- is.na(p) | (p < 0) | (p > 1)
      zero <- (p == 0)
      one  <- (p == 1)
      q    <- numeric(length(p))
      aux1 <- which(na | one | zero)
      if(length(p)!=length(aux1)) {
        aux2 <- c(1:length(p))
        #aux2<- ifelse(length(aux1)!=0,aux2[-c(aux1)],aux2)
        for (i in aux2){
          if (shape==0) q[i]  <- qt(p[i],nu)
          else { qti=optimize(function(x) abs(pstn(x,0,1,shape,nu)-p[i]),c((1.125 - 2.575*sign(shape))*sign(shape),(1.125 + 2.575*sign(shape))*sign(shape)))$minimum
          q[i] <- qti}
        }
      }
    }
    q <- replace(q, na, NA)
    q <- replace(q, zero, -Inf)
    q <- replace(q, one, Inf)
    return(location+sqrt(scale)*q)
  }



qssl <- function (p,location=0, scale=1, shape=0, nu=30, dp=NULL) {
    if(!is.null(dp)) {
      if(!missing(shape))
        stop("You cannot set both component parameters and dp")
      location <- dp[1]
      scale <- dp[2]
      shape <- dp[3]
      nu <- dp[4]
    }
    else
    {
      if (nu<=0) {
        stop("Parameter nu must be positive")
      }
      if (scale<=0) {
        stop("Parameter scale must be positive")
      }
      na <- is.na(p) | (p < 0) | (p > 1)
      zero <- (p == 0)
      one <- (p == 1)
      q <- numeric(length(p))
      aux1 <- which(na | one | zero)
      if(length(p)!=length(aux1)) {
        aux2 <- c(1:length(p))
        # aux2<- ifelse(length(aux1)!=0,aux2[-c(aux1)],aux2)
        for (i in aux2){
          q[i]  <-
            if (shape==0) optimize(function(x) abs(pssl(x,0,1,shape,nu)-p[i]),c(-3,3))$minimum
          else         optimize(function(x) abs(pssl(x,0,1,shape,nu)-p[i]),c((1.125 - 2.575*sign(shape))*sign(shape),(1.125 + 2.575*sign(shape))*sign(shape)))$minimum
        }
      }
    }
    q <- replace(q, na, NA)
    q <- replace(q, zero, -Inf)
    q <- replace(q, one, Inf)
    return(location+sqrt(scale)*q)
  }




qscn <-function (p,location=0, scale=1, shape=0, nu=1, gama=1, dp=NULL) {
    if(!is.null(dp)) {
      if(!missing(shape))
        stop("You cannot set both component parameters and dp")
      location <- dp[1]
      scale    <- dp[2]
      shape    <- dp[3]
      nu       <- dp[4]
      gama     <- dp[5]
    }
    else
    {
      if (nu<=0|nu>=1) {
        stop("Parameter nu must be between 0 and 1.0")
      }
      if(gama<=0|gama>=1) {
        stop("Parameter gama must be between 0 and 1.0")
      }
      if (scale<=0) {
        stop("Parameter scale must be positive")
      }
      na   <- is.na(p) | (p < 0) | (p > 1)
      zero <- (p == 0)
      one  <- (p == 1)
      q    <- numeric(length(p))
      aux1 <- which(na | one | zero)
      if(length(p)!=length(aux1)) {
        aux2 <- c(1:length(p))
        # aux2<- ifelse(length(aux1)!=0,aux2[-c(aux1)],aux2)
        for (i in aux2){
          q[i]  <-
            if (shape==0) optimize(function(x) abs(pscn(x,0,1,shape,nu,gama)-p[i]),c(-3,3))$minimum
          else         optimize(function(x) abs(pscn(x,0,1,shape,nu,gama)-p[i]),c((1.125 - 2.575*sign(shape))*sign(shape),(1.125 + 2.575*sign(shape))*sign(shape)))$minimum
        }
      }
    }
    q <- replace(q, na, NA)
    q <- replace(q, zero, -Inf)
    q <- replace(q, one, Inf)
    return(location+sqrt(scale)*q)
  }



qsep <-function (p,location=0, scale=1, shape=0, nu=1, dp=NULL) {
    if(!is.null(dp)) {
      if(!missing(shape))
        stop("You cannot set both component parameters and dp")
      location <- dp[1]
      scale <- dp[2]
      shape <- dp[3]
      nu <- dp[4]
    }
    else
    {
      if (nu<=0.5|nu>1) {
        stop("Parameter nu must be 0.5 < nu <= 1.0")
      }
      if (scale<=0) {
        stop("Parameter scale must be positive")
      }
      na   <- is.na(p) | (p < 0) | (p > 1)
      zero <- (p == 0)
      one  <- (p == 1)
      q    <- numeric(length(p))
      aux1 <- which(na | one | zero)
      if(length(p)!=length(aux1)) {
        aux2 <- c(1:length(p))
        # aux2<- ifelse(length(aux1)!=0,aux2[-c(aux1)],aux2)
        for (i in aux2){
          q[i]  <-
            if (shape==0) optimize(function(x) abs(psep(x,0,1,shape,nu)-p[i]),c(-3,3))$minimum
          else optimize(function(x) abs(psep(x,0,1,shape,nu)-p[i]),c((1.125 - 2.575*sign(shape))*sign(shape),(1.125 + 2.575*sign(shape))*sign(shape)))$minimum
        }
      }
    }
    q <- replace(q, na, NA)
    q <- replace(q, zero, -Inf)
    q <- replace(q, one, Inf)
    return(location+sqrt(scale)*q)
  }
}

## Auxiliares

psnn<- function (x, location=0, scale=1, shape=0, dp=NULL)
{
  if(!is.null(dp)) {
    if(!missing(shape))
      stop("You cannot set both component parameters and dp")
    location <- dp[1]
    scale <- dp[2]
    shape <- dp[3]
  }
  else
  {
    if(scale<0) {
      stop("Parameter scale must be positive")
    }
    #   z <- (x-location)/sqrt(scale)
    p <- numeric(length(x))
    for (i in 1:length(x)){
      p[i]  <-
        if(abs(x[i]) == Inf)   (1+sign(x[i]))/2
      else{
        integrate(dsnn, -Inf, x[i], location=location, scale=scale, shape = shape)$value
      }
    }
    pmax(0,pmin(1,p))
  }
}

dsnn <- function (x, location=0, scale=1, shape=0, dp=NULL)
{
  if(!is.null(dp)) {
    if(!missing(shape))
      stop("You cannot set both component parameters and dp")
    location <- dp[1]
    scale <- dp[2]
    shape <- dp[3]
  }
  if(scale<0) {
    stop("Parameter scale must be positive")
  }
  if (any(is.na(x))) x<- x[-which(is.na(x))]
  z<-(x-location)/sqrt(scale)
  y<-2*dnorm(z)*pnorm(shape*z)/sqrt(scale)
  return(y)
}


pstn<- function (x, location=0, scale=1, shape=0, nu=30, dp=NULL)
{
  if(!is.null(dp)) {
    if(!missing(shape))
      stop("You cannot set both component parameters and dp")
    location <- dp[1]
    scale <- dp[2]
    shape <- dp[3]
    nu <- dp[4]
  }
  else
  {
    if (nu<0) {
      stop("Parameter nu must be positive")
    }
    if (scale<0) {
      stop("Parameter scale must be positive")
    }
    #   z <- (x-location)/sqrt(scale)
    p <- numeric(length(x))
    for (i in 1:length(x)){
      p[i]  <-
        if(abs(x[i]) == Inf)   (1+sign(x[i]))/2
      else{
        integrate(dstn, -Inf, x[i], location=location, scale=scale, shape = shape, nu = nu)$value
      }
    }
    pmax(0,pmin(1,p))
  }
}

dstn <- function(x, location=0, scale=1, shape=0, nu=30, dp=NULL)
{
  if(!is.null(dp)) {
    if(!missing(shape))
      stop("You cannot set both component parameters and dp")
    location <- dp[1]
    scale <- dp[2]
    shape <- dp[3]
    nu <- dp[4]
  }
  if (nu<0) {
    stop("Parameter nu must be positive")
  }
  if(scale<0) {
    stop("Parameter scale must be positive")
  }
  if (any(is.na(x))) x<- x[-which(is.na(x))]
  z<-(x-location)/sqrt(scale)
  cte<-gamma((nu+1)/2)/gamma(nu/2)/sqrt(nu*pi)
  yt<-cte*(1+z^2/nu)^(-(nu+1)/2)
  y<-2*yt*pnorm(shape*z)/sqrt(scale)
  return(y)
}


pssl<- function(x, location=0, scale=1, shape=0, nu=30, dp=NULL)
{
  if(!is.null(dp)) {
    if(!missing(shape))
      stop("You cannot set both component parameters and dp")
    location <- dp[1]
    scale <- dp[2]
    shape <- dp[3]
    nu<- dp[4]
  }
  else
  {
    if (nu<0) {
      stop("Parameter nu must be positive")
    }
    if(scale<0) {
      stop("Parameter scale must be positive")
    }
    #   z <- (x-location)/sqrt(scale)
    p <- numeric(length(x))
    for (i in 1:length(x)){
      p[i]  <-
        if(abs(x[i]) == Inf)   (1+sign(x[i]))/2
      else{
        integrate(dssl, -Inf, x[i], location=location, scale=scale, shape = shape, nu = nu)$value
      }
    }
    pmax(0,pmin(1,p))
  }
}


dssl <- function(x, location=0, scale=1, shape=0, nu=30, dp=NULL)
{
  if(!is.null(dp)) {
    if(!missing(shape))
      stop("You cannot set both component parameters and dp")
    location <- dp[1]
    scale <- dp[2]
    shape <- dp[3]
    nu <- dp[4]
  }
  if (nu<0) {
    stop("Parameter nu must be positive")
  }
  if(scale<0) {
    stop("Parameter scale must be positive")
  }
  if (any(is.na(x))) x<- x[-which(is.na(x))]
  z<-(x-location)/sqrt(scale)
  n<-length(x)
  cte=nu/((2*pi*scale)^(1/2))
  d=z^2
  fdps=cte*gamma(nu+1/2)*pgamma(1,nu+1/2,scale=2/d)/((d/2)^(nu+1/2))
  fsl=2*fdps*pnorm(shape*z)
  return(fsl)
}


pscn <- function(x, location=0, scale=1, shape=0, nu=1, gama=1, dp=NULL)
{
  if(!is.null(dp)) {
    if(!missing(shape)) {
      stop("You cannot set both component parameters and dp") }
    location <- dp[1]
    scale <- dp[2]
    shape <- dp[3]
    nu <- dp[4]
    gama <- dp[5]
  }
  else
  {
    if (nu<0|nu>1) {
      stop("Parameter nu must be between 0.5 and 1.0")
    }
    if(gama<0|gama>1) {
      stop("Parameter gama must be between 0.5 and 1.0")
    }
    if(scale<0) {
      stop("Parameter scale must be positive")
    }
    #   z <- (x-location)/sqrt(scale)
    p <- numeric(length(x))
    for (i in 1:length(x)){
      p[i]  <-
        if(abs(x[i]) == Inf)   (1+sign(x[i]))/2
      else{
        integrate(dscn, -Inf, x[i], location=location, scale=scale, shape = shape, nu = nu, gama=gama)$value
      }
    }
    pmax(0,pmin(1,p))
  }
}

dscn <- function(x, location=0, scale=1, shape=0, nu=1,gama=1,dp=NULL)
{
  if(!is.null(dp)) {
    if(!missing(shape))
      stop("You cannot set both component parameters and dp")
    location <- dp[1]
    scale <- dp[2]
    shape <- dp[3]
    nu <- dp[4]
    gama <- dp[5]
  }
  if (nu<0|nu>1) {
    stop("Parameter nu must be between 0.5 and 1.0")
  }
  if(gama<0|gama>1) {
    stop("Parameter gama must be between 0.5 and 1.0")
  }
  if(scale<0) {
    stop("Parameter scale must be positive")
  }
  if (any(is.na(x))) x<- x[-which(is.na(x))]
  z<-(x-location)/sqrt(scale)
  z2<-z*sqrt(gama)
  y<-2*(nu*dnorm(x,mean=location,sd=sqrt(scale/gama))+(1-nu)*dnorm(x,mean=location,sd=sqrt(scale)))*pnorm(shape*z)
  return(y)
}


psep <- function(x, location=0, scale=1, shape=0, nu=1, dp=NULL)
{
  if(!is.null(dp)) {
    if(!missing(shape)) {
      stop("You cannot set both component parameters and dp") }
    location <- dp[1]
    scale <- dp[2]
    shape <- dp[3]
    nu <- dp[4]
  }
  else
  {
    if (nu<0.5|nu>1) {
      stop("Parameter nu must be between 0.5 and 1.0")
    }
    if(scale<0) {
      stop("Parameter scale must be positive")
    }
    #   z <- (x-location)/sqrt(scale)
    p <- numeric(length(x))
    for (i in 1:length(x)){
      p[i]  <-
        if(abs(x[i]) == Inf)   (1+sign(x[i]))/2
      else{
        integrate(dsep, -Inf, x[i], location=location, scale=scale, shape = shape, nu = nu)$value
      }
    }
    pmax(0,pmin(1,p))
  }
}

dsep <- function(x, location=0, scale=1, shape=0, nu=1, dp=NULL)
{
  if(!is.null(dp)) {
    if(!missing(shape))
      stop("You cannot set both component parameters and dp")
    location <- dp[1]
    scale <- dp[2]
    shape <- dp[3]
    nu<- dp[4]
  }
  if (nu<0.5|nu>1) {
    stop("Parameter nu must be between 0.5 and 1.0")
  }
  if(scale<0) {
    stop("Parameter scale must be positive")
  }
  if (any(is.na(x))) x<- x[-which(is.na(x))]
  nu1<-1/(2*nu)
  z<-(x-location)/sqrt(scale)
  d<-z^2
  cte<-nu/(2^nu1)/gamma(nu1)/sqrt(scale)
  yt<-cte*exp(-0.5*d^(nu))
  y<-2*yt*pnorm(shape*z)
  return(y)
}


################################################################################
# Random Generation of SSMN distributions

rssmn <- function(n, location=0, scale=1, shape=0, nu= 1, gama=1, dp=NULL, family="sn")
{
 if(family == "sn")
 {
  if(!is.null(dp))
  {
   if(!missing(shape)) {stop("You cannot set both component parameters and dp")}
   location <- dp[1]
   scale    <- dp[2]
   shape    <- dp[3]
  }
  if(scale<=0) {stop("Parameter scale must be positive")}
  delta <- shape/sqrt(1+shape^2)
  u1    <- rnorm(n)
  u2    <- rnorm(n)
  y     <- location+sqrt(scale)*(delta*abs(u1)+sqrt(1-delta^2)*u2)
 }

 if(family == "stn")
 {
   #rstn <- function(n, location, scale, shape, nu, dp=NULL)
   #{
    if(!is.null(dp))
    {
     if(!missing(shape)){stop("You cannot set both component parameters and dp")}
     location <- dp[1]
     scale    <- dp[2]
     shape    <- dp[3]
     nu       <- dp[4]
    }
    if (nu<=0) {stop("Parameter nu must be positive")}
    if(scale<=0){stop("Parameter scale must be positive")}
    if (nu<1){warning('Nu < 1 can generate values tending to infinite',call. = FALSE)}
    u      <- rgamma(n,nu/2,nu/2)
    ku     <- 1/u
    shape1 <- shape*sqrt(ku)
    scale1 <- scale*ku
    delta  <- shape1/sqrt(1+shape1^2)
    u1     <- rnorm(n)
    u2     <- rnorm(n)
    y      <- location+sqrt(scale1)*(delta*abs(u1)+sqrt(1-delta^2)*u2)
   #return(y)
   #}

  }

 if(family == "ssl")
 {
  #rssl <- function(n, location, scale, shape, nu, dp)
  #{
   if(!is.null(dp))
   {
    if(!missing(shape)){stop("You cannot set both component parameters and dp")}
    location <- dp[1]
    scale    <- dp[2]
    shape    <- dp[3]
    nu       <- dp[4]
   }
   if(nu<=0){stop("Parameter nu must be positive")}
   if(scale<0){stop("Parameter scale must be positive")}
   v      <- runif(n,0,1)
   u      <- v^(1/nu)
   ku     <- 1/u
   shape1 <- shape*sqrt(ku)
   scale1 <- scale*ku
   delta  <- shape1/sqrt(1+shape1^2)
   u1     <- rnorm(n)
   u2     <- rnorm(n)
   y      <- location+sqrt(scale1)*(delta*abs(u1)+sqrt(1-delta^2)*u2)
   #return(y)
   #}
 }

 if(family == "scn")
 {
   #rscn <- function(n, location, scale, shape, nu, gama, dp)
   #{
     if(!is.null(dp)) {
       if(!missing(shape))
         stop("You cannot set both component parameters and dp")
       location <- dp[1]
       scale <- dp[2]
       shape <- dp[3]
       nu <- dp[4]
       gama <- dp[5]
     }
     if (nu<0|nu>1) {
       stop("Parameter nu must be between 0 and 1.0")
     }
     if(gama<=0|gama>1) {
       stop("Parameter gama must be between 0 and 1.0")
     }
     if(scale<=0) {
       stop("Parameter scale must be positive")
     }
     uu     <- rbinom(n,1,nu)
     u      <- gama*uu+1-uu
     ku     <- 1/u
     shape1 <- shape*sqrt(ku)
     scale1 <- scale*ku
     delta  <- shape1/sqrt(1+shape1^2)
     u1     <- rnorm(n)
     u2     <- rnorm(n)
     y      <- location+sqrt(scale1)*(delta*abs(u1)+sqrt(1-delta^2)*u2)
    # return(y)
   #}
 }

 if(family == "sep")
 {
  #rsep <- function (n, location, scale, shape, nu, dp)
  #{
   if(!is.null(dp))
   {
    if(!missing(shape)) {stop("You cannot set both component parameters and dp")}
     location <- dp[1]
     scale    <- dp[2]
     shape    <- dp[3]
     nu       <- dp[4]
     }
     if (nu<=0.5|nu>1){stop("Parameter nu must be 0.5 < nu <= 1.0")}
     if(scale<=0) {stop("Parameter scale must be positive")}
    t1  <- rgamma(n,shape=nu/2,scale=2)
    R   <- t1^(1/(2*nu))
    u   <- sign(rnorm(n))
    yep <- u*R
    w   <- rnorm(n)
    y   <- location+sqrt(scale)*yep*sign(shape*yep-w)
    #return(y)
   #}
 }

 return(list(y=y,family=family))

}

#rssmn(n=20, location=0, scale=1, shape=0, nu= 1, gama=1, dp=NULL, family="sn")


###############################################################################
# Negative of Log-Likelihood of SSMN Regression Models

ssmn.logL <- function(theta,y,X, family="sn")
{
  if(family=="sn")  {
    # Theta=[beta,sigma2,lambda]
    if(missing(X)) {X<-matrix(1,n,1)}
    p<-signif(ncol(X),digits=7)
    beta<-signif(theta[1:p],digits=7)
    sigma2<-signif(theta[1+p],digits=7)
    sigma<-signif(sqrt(sigma2),digits=7)
    lambda<-signif(theta[p+2],digits=7)
    vero<-signif(2*dnorm(y,X%*%beta,sigma)*pnorm(lambda*(y-X%*%beta)/sigma),digits=7)
    #logvero<-signif(log(vero),digits=7)
    logvero=-sum(log(vero))
  }
  if(family=="stn") {
    # Theta=[beta,sigma2,lambda,nu]
    if(missing(X)) {X<-matrix(1,n,1)}
    p<-ncol(X)
    beta<-signif(theta[1:p],digits=7)
    sigma2<-signif(theta[p+1],digits=7)
    lambda<-signif(theta[p+2],digits=7)
    nu<-signif(theta[p+3],digits=7)
    z=(y-X%*%beta)/sqrt(sigma2)
    d<-z^2
    aux<-signif(lambda*z,digits=7)
    aux1<-signif(pmax(aux,-37),digits=7)
    cnugama<-signif(2*gamma((nu+1)/2)/(gamma(nu/2)*sqrt(nu*pi)),digits=7)
    vero<-signif(cnugama/sqrt(sigma2)*(1+d/nu)^(-(nu+1)/2)*pnorm(aux1),digits=7)
    logvero=-sum(log(vero))
  }
  if(family=="ssl") {
    #% =Theta=[beta,sigma^2,lambda,nu]
    n<-length(y)
    if(missing(X)) {X<-matrix(1,n,1)}
    p<-ncol(X)
    beta<-theta[1:p]
    sigma2<-theta[p+1]
    lambda<-theta[p+2]
    nu<-theta[p+3]
    z=(y-X%*%beta)/sqrt(sigma2)
    d<-z^2
    aux<-signif(lambda*z,digits=7)
    aux1<-signif(pmax(aux,-37),digits=7)
    cte=nu/((2*pi*sigma2)^0.5)
    fdps=cte*gamma(nu+0.5)*pgamma(1,nu+1/2,scale=2/d)/((d/2)^(nu+0.5))
    fdps[which(z==0)]=cte/(nu+0.5)
    vero=2*fdps*pnorm(aux1)
    logvero=-sum(log(vero))
  }
  if(family=="scn") {
    #% Theta=[beta,sigma2,lambda,nu,gama]
    if(missing(X)) {X<-matrix(1,n,1)}
    p<-ncol(X)
    beta<-theta[1:p]
    mu=X%*%beta
    sigma2<-theta[p+1]
    lambda<-theta[p+2]
    nu<-theta[p+3]
    gama<-theta[p+4]
    z=(y-X%*%beta)/sqrt(sigma2)
    aux<-signif(lambda*z,digits=7)
    aux1<-signif(pmax(aux,-37),digits=7)
    vero<-signif(2*(nu*dnorm(y,mu,sqrt(sigma2/gama))+(1-nu)*dnorm(y,mu,sqrt(sigma2)))*pnorm(aux1),digits=7)
    logvero=-sum(log(vero))
  }
  if(family=="sep") {
    #% Theta=[beta,sigma2,lambda,nu]
    p<-ncol(X)
    beta<-theta[1:p]
    sigma2<-theta[p+1]
    lambda<-theta[p+2]
    nu<-theta[p+3]
    z=(y-X%*%beta)/sqrt(sigma2)
    d<-z^2
    aux<-signif(lambda*z,digits=7)
    aux1<-signif(pmax(aux,-37),digits=7)
    cnu<-2*nu/(2^(1/(2*nu))*sqrt(sigma2)*gamma(1/(2*nu)))
    vero<-cnu*exp(-0.5*d^nu)*pnorm(aux1)
    logvero=-sum(log(vero))
  }
  return(logvero)
}


############################################################################
# Observed Fisher's Information Matrix

ssmn.info <- function(y,x,theta, family="sn")
{
  if(family=="sn")  {
    # THETA: (beta,sigma2,lambda)
    n=length(y)
    if(missing(x)) {x<-matrix(1,n,1)}
    p<-ncol(x)
    beta<-theta[1:p]
    sigma2<-as.numeric(theta[p+1])
    lambda<-as.numeric(theta[p+2])
    res<-y-x%*%beta
    res<-as.vector(res)
    sigma<-sqrt(sigma2)
    aux<-lambda*res/sigma
    aux1<-pmax(aux,-37)
    Wphi<-dnorm(aux1)/pnorm(aux1)
    d<-res^2/sigma2
    ######################### I2(beta,sigma2,lambda) ###########################
    Qbeta<-t(res)%*%res
    Wphi1<- -Wphi*(aux1+Wphi)
    Qwbeta<- t(res)%*%diag(Wphi1)%*%res
    I2<-matrix(0,p+2,p+2)
    I2[1:p,1:p]<- -lambda^2/sigma2*t(x)%*%diag(Wphi1)%*%x
    I2[1:p,p+1]<- -lambda/(2*sigma^3)*t(x)%*%Wphi-lambda^2/(2*sigma2^2)*t(x)%*%diag(Wphi1)%*%res
    I2[1:p,p+2]<- 1/sigma*t(x)%*%Wphi+lambda/sigma2*t(x)%*%diag(Wphi1)%*%res
    I2[p+1,p+1]<- -3*lambda/(4*sigma^5)*t(res)%*%Wphi-lambda^2/(4*sigma2^3)*Qwbeta
    I2[p+1,p+2]<- 1/(2*sigma^3)*t(res)%*%Wphi+lambda/(2*sigma2^2)*Qwbeta
    I2[p+2,p+2]<- -Qwbeta/sigma2
    I2[p+2,p+1]<- I2[p+1,p+2]
    I2[p+1,1:p]<- t(I2[1:p,p+1])
    I2[p+2,1:p]<- t(I2[1:p,p+2])
    ######################## I1(beta,sigma2,lambda) ###########################
    I1<-matrix(0,p+2,p+2)
    I1[1:p,1:p]<-1/(sigma2)*t(x)%*%x
    I1[1:p,p+1]<-1/(sigma2^2)*t(x)%*%res
    I1[p+1,p+1]<- -n/(2*sigma2^2)+1/(sigma2^3)*t(res)%*%res
    I1[p+1,1:p]<-t(I1[1:p,p+1])
    ########
    Info <- I1+I2
  }
  if(family=="stn") {
    # THETA: (beta,sigma2,lambda,nu)
    n<-nrow(x)
    if(missing(x)) {x<-matrix(1,n,1)}
    p<-ncol(x)
    beta<-theta[1:p]
    sigma2<-as.numeric(theta[p+1])
    lambda<-as.numeric(theta[p+2])
    nu<-as.numeric(theta[p+3])
    mu<-x%*%beta
    res<-as.vector(y-mu)
    sigma<-sqrt(sigma2)
    aux<-lambda*res/sigma
    aux1<-pmax(aux,-37)
    Wphi<-dnorm(aux1)/pnorm(aux1)
    d<-res^2/sigma2
    ######################### I2(beta,sigma2,lambda) ###########################
    Qbeta<-t(res)%*%res
    Wphi1<- as.vector(-Wphi*(aux1+Wphi))
    Qwbeta<-t(res)%*%diag(Wphi1)%*%res
    I2<-matrix(0,p+3,p+3)
    I2[1:p,1:p]<- -lambda^2/sigma2*t(x)%*%diag(Wphi1)%*%x
    I2[1:p,p+1]<- -lambda/(2*sigma^3)*t(x)%*%Wphi-lambda^2/(2*sigma2^2)*t(x)%*%diag(Wphi1)%*%res
    I2[1:p,p+2]<- 1/sigma*t(x)%*%Wphi+lambda/sigma2*t(x)%*%diag(Wphi1)%*%res
    I2[p+1,p+1]<- -3*lambda/(4*sigma^5)*t(res)%*%Wphi-lambda^2/(4*sigma2^3)*Qwbeta
    I2[p+1,p+2]<- 1/(2*sigma^3)*t(res)%*%Wphi+lambda/(2*sigma2^2)*Qwbeta
    I2[p+2,p+2]<- -Qwbeta/sigma2
    I2[p+2,p+1]<- I2[p+1,p+2]
    I2[p+1,1:p]<- t(I2[1:p,p+1])
    I2[p+2,1:p]<- t(I2[1:p,p+2])
    ######################## I1(beta,sigma2,lambda,nu) ###########################
    I1<-matrix(0,p+3,p+3)
    V<-(1+d/nu)^(-1)
    Qvbeta<- t(res)%*%diag(V)%*%res
    hnu<-psigamma((nu+1)/2,1)-psigamma(nu/2,1)+2/(nu^2)
    I1[1:p,1:p]<- -(nu+1)/(nu*sigma2)*t(x)%*%diag(V)%*%(2/nu*diag(V)%*%diag(d)-diag(n))%*%x
    I1[1:p,p+1]<- -(nu+1)/(nu*sigma2^2)*t(x)%*%diag(res)%*%(1/nu*diag(V)%*%diag(d)-diag(n))%*%V
    I1[1:p,p+3]<- -1/(nu^2*sigma2)*t(x)%*%diag(res)%*%((nu+1)/nu*diag(V)%*%diag(d)-diag(n))%*%V
    I1[p+1,p+1]<- -n/(2*(sigma2^2))+(nu+1)/(nu*sigma2^3)*Qvbeta-(nu+1)/(2*(nu^2)*(sigma2^4))*t(res)%*%(diag(res^3))%*%diag(V)%*%V
    I1[p+1,p+3]<- -0.5/((nu*sigma2)^2)*t(res)%*%diag(V)%*%((nu+1)/nu*diag(V)%*%diag(d)-diag(n))%*%res
    I1[p+3,p+3]<- -n/4*hnu-0.5/(nu^2*sigma2)*Qvbeta-0.5/((nu^4)*sigma2)*t(res)%*%diag(V)%*%((nu+1)*diag(V)%*%diag(d)-nu*(nu+2)*diag(n))%*%res
    I1[p+3,1:p]<-I1[1:p,p+3]
    I1[p+1,1:p]<-I1[1:p,p+1]
    I1[p+3,p+1]<-I1[p+1,p+3]
    ########
    I<-I1+I2
    return(I)
  }
  if(family=="ssl") {
    # THETA: (beta,sigma2,lambda,nu)
    n<-nrow(x)
    if(missing(x)) {x<-matrix(1,n,1)}
    p<-ncol(x)
    beta<-theta[1:p]
    sigma2<-as.numeric(theta[p+1])
    lambda<-as.numeric(theta[p+2])
    nu<-as.numeric(theta[p+3])
    res<-as.vector(y-x%*%beta)
    sigma<-sqrt(sigma2)
    aux<-lambda*res/sigma
    aux1<-pmax(aux,-37)
    Wphi<-dnorm(aux1)/pnorm(aux1)
    d<-res^2/sigma2
    ######################### I2(beta,sigma2,lambda) ###########################
    Qbeta<-t(res)%*%res
    Wphi1<- as.vector(-Wphi*(aux1+Wphi))
    Qwbeta<-t(res)%*%diag(Wphi1)%*%res
    I2<-matrix(0,p+3,p+3)
    I2[1:p,1:p]<- -lambda^2/sigma2*t(x)%*%diag(Wphi1)%*%x
    I2[1:p,p+1]<- -lambda/(2*sigma^3)*t(x)%*%Wphi-lambda^2/(2*sigma2^2)*t(x)%*%diag(Wphi1)%*%res
    I2[1:p,p+2]<- 1/sigma*t(x)%*%Wphi+lambda/sigma2*t(x)%*%diag(Wphi1)%*%res
    I2[p+1,p+1]<- -3*lambda/(4*sigma^5)*t(res)%*%Wphi-lambda^2/(4*sigma2^3)*Qwbeta
    I2[p+1,p+2]<- 1/(2*sigma^3)*t(res)%*%Wphi+lambda/(2*sigma2^2)*Qwbeta
    I2[p+2,p+2]<- -Qwbeta/sigma2
    I2[p+2,p+1]<- I2[p+1,p+2]
    I2[p+1,1:p]<- t(I2[1:p,p+1])
    I2[p+2,1:p]<- t(I2[1:p,p+2])
    ######################## I1(beta,sigma2,lambda,nu) ###########################
    I1<-matrix(0,p+3,p+3)
 P11<-pgamma(1,nu+0.5,scale=2/d)
 P13<-pgamma(1,nu+1.5,scale=2/d)
 P15<-pgamma(1,nu+2.5,scale=2/d)
 IG1<-gamma(nu+0.5)*P11/((d/2)^(nu+0.5))
 IG3<-gamma(nu+1.5)*P13/((d/2)^(nu+1.5))
 IG5<-gamma(nu+2.5)*P15/((d/2)^(nu+2.5))
 rp<-5000
    u1<-matrix(0,n,rp)
    u2<-u1
    for (j in 1:n) {
       u1[j,]<-rtrunc(rp,"gamma",a=0,b=1,shape=nu+1/2,rate=d[j]/2)
       u2[j,]<-rtrunc(rp,"gamma",a=0,b=1,shape=nu+3/2,rate=d[j]/2)
    }
    #E11=gamma(nu+0.5)*((d/2)^(-(nu+0.5)))*P11*rowMeans(log(u1))
    #E13=gamma(nu+1.5)*((d/2)^(-(nu+1.5)))*P13*rowMeans(log(u2))
    #E21=gamma(nu+2.5)*((d/2)^(-(nu+2.5)))*P11*rowMeans(log(u1)^2)
    E11=P11*rowMeans(log(u1))
    E13=P13*rowMeans(log(u2))
    E21=P11*rowMeans(log(u1)^2)
    I1[1:p,1:p]<- -1/sigma2*t(x)%*%diag(-IG3/IG1+d*(IG5-IG3^2/IG1)/IG1)%*%x
    I1[1:p,p+1]<- -1/(2*sigma2^2)*t(x)%*%diag(res/IG1)%*%(-2*IG3+d*(IG5-IG3^2/IG1))
    I1[1:p,p+3]<- -1/(sigma2)*t(x)%*%diag(res/(IG1^2))%*%(IG1*E13-IG3*E11)
    I1[p+1,p+1]<- -1/(2*sigma2^2)*t(n+(d/IG1))%*%(-2*IG3+0.5*d*(IG5-IG3^2/IG1))
    I1[p+1,p+3]<- -1/(2*sigma2)*t(d/(IG1^2))%*%(IG1*E13-IG3*E11)
    I1[p+3,p+3]<- n/(nu^2)-t(IG1*E21-E11^2)%*%(IG1^(-2))
    I1[p+1,1:p]<- t(I1[1:p,p+1])
    I1[p+3,1:p]<- t(I1[1:p,p+3])
    I1[p+3,p+1]<- I1[p+1,p+3]
    ########
    I<-I1+I2
    return(I)
  }
  if(family=="scn") {
    # THETA: (beta,sigma2,lambda,nu,gama)
    n<-nrow(x)
    if(missing(x)) {x<-matrix(1,n,1)}
    p<-ncol(x)
    beta<-theta[1:p]
    sigma2<-as.numeric(theta[p+1])
    lambda<-as.numeric(theta[p+2])
    nu<-as.numeric(theta[p+3])
    gama<-as.numeric(theta[p+4])
    res<-y-x%*%beta
    sigma<-sqrt(sigma2)
    aux<-lambda*res/sigma
    aux1<-pmax(aux,-37)
    Wphi<-dnorm(aux1)/pnorm(aux1)
    d<-res^2/sigma2
    ######################### I2(beta,sigma2,lambda) ###########################
    Qbeta<-t(res)%*%res
    Wphi1<- as.vector(-Wphi*(aux1+Wphi))
    Qwbeta<-t(res)%*%diag(Wphi1)%*%res
    I2<-matrix(0,p+4,p+4)
    I2[1:p,1:p]<- -lambda^2/sigma2*t(x)%*%diag(Wphi1)%*%x
    I2[1:p,p+1]<- -lambda/(2*sigma^3)*t(x)%*%Wphi-lambda^2/(2*sigma2^2)*t(x)%*%diag(Wphi1)%*%res
    I2[1:p,p+2]<- 1/sigma*t(x)%*%Wphi+lambda/sigma2*t(x)%*%diag(Wphi1)%*%res
    I2[p+1,p+1]<- -3*lambda/(4*sigma^5)*t(res)%*%Wphi-lambda^2/(4*sigma2^3)*Qwbeta
    I2[p+1,p+2]<- 1/(2*sigma^3)*t(res)%*%Wphi+lambda/(2*sigma2^2)*Qwbeta
    I2[p+2,p+2]<- -Qwbeta/sigma2
    I2[p+2,p+1]<- I2[p+1,p+2]
    I2[p+1,1:p]<- t(I2[1:p,p+1])
    I2[p+2,1:p]<- t(I2[1:p,p+2])
    ######################## I1(beta,sigma2,lambda,nu,gama) #######################
    I1<-matrix(0,p+4,p+4)
    sigma<-sqrt(sigma2)
    fi<-dnorm(y,x%*%beta,sigma)
    figama<-dnorm(y,x%*%beta,sqrt(sigma2/gama))
    fs<-nu*figama+(1-nu)*fi
    aux2<-as.vector(res*(nu*gama*figama+(1-nu)*fi))
    dfsbeta<-1/sigma2*t(diag(aux2)%*%x)
    dlbeta<-matrix(0,p,n)
    for (i in 1:n) {
      dlbeta[,i]<-dfsbeta[,i]/fs[i]
    }
    dfssigma<-1/(2*sigma2)*(nu*figama*(gama*d-1)+(1-nu)*fi*(d-1))
    dfsnu<-figama-fi
    dfsgama<-nu/(2*gama)*figama*(1-gama*d)
    d2lbetabeta<-matrix(0,p,p)
    for (i in 1:n) {
      auxn<-nu*gama*figama[i]+(1-nu)*fi[i]
      aux<-1/(sigma2*fs[i])*(d[i]*(nu*(gama^2)*figama[i]+(1-nu)*fi[i])-auxn-d[i]/fs[i]*auxn^2)*(as.vector(x[i,]))%*%t(as.vector(x[i,]))
      d2lbetabeta<-d2lbetabeta+aux
    }
    d2lsigmabeta<-matrix(0,p,1)
    for (i in 1:n) {
      d2sigmabeta<-res[i]/(2*sigma2^2)*(nu*gama*figama[i]*(gama*d[i]-3)+(1-nu)*fi[i]*(d[i]-3))*as.vector(x[i,])
      aux<-(d2sigmabeta*fs[i]-dfssigma[i]*dfsbeta[,i])/(fs[i]^2)
      d2lsigmabeta<-d2lsigmabeta+aux
    }
    d2lnubeta<-matrix(0,p,1)
    for (i in 1:n) {
      d2nubeta<-1/sigma2*(gama*figama[i]-fi[i])*res[i]*as.vector(x[i,])
      aux<--1/(fs[i]^2)*dfsnu[i]*dfsbeta[,i]+1/fs[i]*d2nubeta
      d2lnubeta<-d2lnubeta+aux
    }
    d2lgamabeta<-matrix(0,p,1)
    for (i in 1:n) {
      d2gamabeta<-nu/(2*sigma2)*figama[i]*(3-gama*d[i])*res[i]*as.vector(x[i,])
      aux<--1/(fs[i]^2)*dfsgama[i]*dfsbeta[,i]+1/fs[i]*d2gamabeta
      d2lgamabeta<-d2lgamabeta+aux;
    }
    d2fssigmasigma=1/(2*sigma2^2)*(-2*nu*gama*figama*d-2*(1-nu)*fi*d+fs+nu/2*figama*(gama*d-1)^2+(1-nu)/2*fi*(d-1)^2);
    Isigmasigma<-1/(fs^2)*(d2fssigmasigma*fs-dfssigma^2)
    d2fsnusigma<-1/(2*sigma2)*(d*(gama*figama-fi)-figama+fi)
    Inusigma<-1/(fs^2)*(d2fsnusigma*fs-dfssigma*dfsnu)
    d2fsgamasigma<-nu/(4*sigma2)*figama*(4*d-gama*(d^2)-1/gama)
    Igamasigma<-1/(fs^2)*(d2fsgamasigma*fs-dfssigma*dfsgama)
    d2fsnunu<-0
    Inunu<- -(dfsnu/fs)^2
    d2fsnugama<-1/2*figama*(1/gama-d)
    Inugama<-1/(fs^2)*(d2fsnugama*fs-dfsnu*dfsgama)
    d2fsgamagama<-nu/4*figama*((1/gama-d)^2-2/(gama^2))
    Igamagama<-1/(fs^2)*(d2fsgamagama*fs-dfsgama^2)
    I1[1:p,1:p]<--d2lbetabeta
    I1[1:p,p+1]<--d2lsigmabeta
    I1[1:p,p+3]<--d2lnubeta
    I1[1:p,p+4]<--d2lgamabeta
    I1[p+1,p+1]<--sum(Isigmasigma)
    I1[p+1,p+3]<--sum(Inusigma)
    I1[p+1,p+4]<--sum(Igamasigma)
    I1[p+3,p+3]<--sum(Inunu)
    I1[p+3,p+4]<--sum(Inugama)
    I1[p+4,p+4]<--sum(Igamagama)
    I1[p+1,1:p]<- t(I1[1:p,p+1])
    I1[p+3,1:p]<- t(I1[1:p,p+3])
    I1[p+4,1:p]<- t(I1[1:p,p+4])
    I1[p+3,p+1]<- I1[p+1,p+3]
    I1[p+4,p+1]<- I1[p+1,p+4]
    I1[p+4,p+3]<- I1[p+3,p+4]
    ########
    I<-I1+I2
    return(I)
  }
  if(family=="sep") {
    # THETA: (beta,sigma2,lambda,nu)
    n<-nrow(x)
    if(missing(x)) {x<-matrix(1,n,1)}
    p<-ncol(x)
    beta<-theta[1:p]
    sigma2<-as.numeric(theta[p+1])
    lambda<-as.numeric(theta[p+2])
    nu<-as.numeric(theta[p+3])
    res<-as.vector(y-x%*%beta)
    sigma<-sqrt(sigma2)
    aux<-lambda*res/sigma
    aux1<-pmax(aux,-37)
    Wphi<-dnorm(aux1)/pnorm(aux1)
    d<-res^2/sigma2
    ################### I2(beta,sigma2,lambda) ########################
    Qbeta<-t(res)%*%res
    Wphi1<- as.vector(-Wphi*(aux1+Wphi))
    Qwbeta<-t(res)%*%diag(Wphi1)%*%res
    I2<-matrix(0,p+3,p+3)
    I2[1:p,1:p]<- -lambda^2/sigma2*t(x)%*%diag(Wphi1)%*%x
    I2[1:p,p+1]<- -lambda/(2*sigma^3)*t(x)%*%Wphi-lambda^2/(2*sigma2^2)*t(x)%*%diag(Wphi1)%*%res
    I2[1:p,p+2]<- 1/sigma*t(x)%*%Wphi+lambda/sigma2*t(x)%*%diag(Wphi1)%*%res
    I2[p+1,p+1]<- -3*lambda/(4*sigma^5)*t(res)%*%Wphi-lambda^2/(4*sigma2^3)*Qwbeta
    I2[p+1,p+2]<- 1/(2*sigma^3)*t(res)%*%Wphi+lambda/(2*sigma2^2)*Qwbeta
    I2[p+2,p+2]<- -Qwbeta/sigma2
    I2[p+2,p+1]<- I2[p+1,p+2]
    I2[p+1,1:p]<- t(I2[1:p,p+1])
    I2[p+2,1:p]<- t(I2[1:p,p+2])
    ######################## I1(beta,sigma2,lambda,nu) ###########################
    I1<-matrix(0,p+3,p+3)
    dnu<- d^(nu-1)
    I1[1:p,1:p]<- nu*(2*nu-1)/sigma2*t(x)%*%diag(dnu)%*%x
    I1[1:p,p+1]<- nu^2/(sigma2^2)*t(x)%*%diag(res)%*%(dnu)
    I1[1:p,p+3]<- -1/(sigma2)*t(x)%*%diag(res*dnu)%*%(matrix(1,n,1)+nu*log(d))
    I1[p+1,p+1]<- 1/(2*(sigma2^2))*(nu*(nu+1)*sum(d^nu)-n)
    I1[p+1,p+3]<- -1/(2*sigma2)*t(d^nu)%*%(matrix(1,n,1)+nu*log(d))
    nu1<- 1/(2*nu)
    I1[p+3,p+3]<- n/(nu^2)*(1+(log(2)+psigamma(nu1))/nu+(nu1^2)*psigamma(nu1,1))+1/2*t(d^nu)%*%((log(d))^2)
    I1[p+1,1:p]<- t(I1[1:p,p+1])
    I1[p+3,1:p]<- t(I1[1:p,p+3])
    I1[p+3,p+1]<- I1[p+1,p+3]
    ########
    I<-I1+I2
    return(I)
  }
  return(Info)
}
#ssmn.info(y,x,theta, family="sn")



gamtruncrnd<- function(a,b,t,n,m)
{
  # Referencia: ANNE PHILIPPE. Simulation of right and left truncated gamma
  # distributions by mixtures. Statistics and Computing (1997) 7, 173-181

  # a: parameter of location
  # t: right truncation, f_t varies from 0 to t
  # [n,m]: sample sample
  # If X ~ TG(a,b,t) then X/t ~ TG(a,bt,1)
  # N: number of replications to ensure that the probability of acceptance algorithm ## acceptance-rejection P(N) > p.
  tp=qnorm(0.99,0,1)
  N=ceiling((tp+sqrt(tp^2+4*b))^2/4) # ceiling(X): integer greater than X
  N=min(N,171)
  K=matrix(1:N)
  M=1/sum(b^(K-1)/gamma(K))
  Y=matrix(0,n,m)
  for (j in 1:n)
  {
    for (k in 1:m)
    {
      u=runif(1,0,1)
      x=betamixt(a,b,N)
      aux=((1-x)*b)^(K-1)/gamma(K)
      rhox=1/(exp(b*x)*sum(aux))
      Y[j,k]=x
      while (u*M > rhox)
      {
        u=runif(1,0,1)
        x=betamixt(a,b,N)
        aux=((1-x)*b)^(K-1)/gamma(K)
        rhox=1/(exp(b*x)*sum(aux))
        Y[j,k]=x
      }
    }
  }
  return(Y)
}

betamixt<- function (a,b,N)
{
  # gera de uma distribui?ao de misturas de betas
  wb=matrix(1,N,1)
  wt=wb
  for (j in 2:N)
  {
    wb[j]=wb[j-1]*b/(a+j-1)
  }
  for (j in 2:N)
  {
    wt[j]=wt[j-1]+wb[j]
  }
  wt=rbind(0,wt/wt[N])
  u=runif(1,0,1)
  dif=wt-u
  k=length(which(dif<0))
  Y=rbeta(1,a,k)
  return(Y)
}



################################################################################
# Negative of Log-likelihood of SMN regression models

ftn<-function(nu,y,X,theta) {
  n<-length(y)
  p=ncol(X)
  beta=theta[1:p]
  sigma2=theta[p+1]
  erro=y-X%*%beta
  aux=1+erro^2/(nu*sigma2)
  logfy=n*log(gamma((nu+1)/2))-n/2*log(nu)-n*log(gamma(nu/2))-(nu+1)/2*sum(log(aux))
  return(-logfy)
}

fsl<-function(nu,y,X,theta){
  n<-length(y)
  p=ncol(X)
  beta=theta[1:p]
  sigma2=theta[p+1]
  cte=nu/sqrt(2*pi*sigma2)*gamma(nu+0.5)
  d=(y-X%*%beta)^2/sigma2
  fy=cte*pgamma(1,nu+0.5,scale=2/d)/((d/2)^(nu+0.5))
  fy[which(d==0)]=cte/(nu+1/2)
  return(-sum(log(fy)))
}

fcn      <- function(nugama,y,X,theta)
{
  n      <- length(y)
  p      <- ncol(X)
  beta   <- theta[1:p]
  mu     <- X%*%beta
  sigma2 <- theta[p+1]
  nu     <- nugama[1]
  gama   <- nugama[2]
  fy     <- nu*dnorm(y,mu,sqrt(sigma2/gama))+(1-nu)*dnorm(y,mu,sqrt(sigma2))
  return(-sum(log(fy)))
}

fep      <- function(nu,y,X,theta0)
{
  n      <- length(y)
  p      <- ncol(X)
  beta   <- theta0[1:p]
  sigma2 <- theta0[p+1]
  d      <- (y-X%*%beta)^2/sigma2
  cnu    <- nu/(2^(1/(2*nu))*sqrt(sigma2)*gamma(1/(2*nu)))
  fy     <- cnu*exp(-0.5*d^nu)
  return(-sum(log(fy)))

  envelSCN <- function(Y,mu,Sigma,nu,gama){
    n=length(Y)
    p=1
    d2=(Y-mu)^2/Sigma
    d2s=sort(d2)
    xq2 <- qchisq(ppoints(n), p)
    Xsim<-matrix(0,200,n)
    for(i in 1:200){
      u1<-rchisq(n, p)/gama
      u2<-rchisq(n, p)
      u3<-runif(n)
      id<-(u3<nu)
      u2[id]<-u1[id]
      Xsim[i,]=u2}
    Xsim2<-apply(Xsim,1,sort)
    d21<-rep(0,n)
    d22<-rep(0,n)
    for(i in 1:n){
      d21[i]  <- quantile(Xsim2[i,],0.025)
      d22[i]  <- quantile(Xsim2[i,],0.975)}
    d2med <-apply(Xsim2,1,mean)

    fy <- range(d2s,d21,d22)
    plot(xq2,d2s,xlab = expression(paste("Theoretical ",chi[1]^2, " quantiles")),
         ylab="Sample values and simulated envelope",pch=20,ylim=fy)
    par(new=T)
    plot(xq2,d21,type="l",ylim=fy,xlab="",ylab="")
    par(new=T)
    plot(xq2,d2med,type="l",ylim=fy,xlab="",ylab="",lty="dashed")
    par(new=T)
    plot(xq2,d22,type="l",ylim=fy,xlab="",ylab="")
  }

}


################################################################################
# Negative of log-likelihood of SSMN regression models, estimating nu/gama
logL <- function(theta,y,X, family)
{
 if(family=="sn")
 {
  #% Theta=[beta,sigma2,lambda]
  p       <- ncol(X)
  beta    <- theta[1:p]
  sigma2  <- theta[p+1]
  lambda  <- theta[p+2]
  z       <- (y-X%*%beta)/sqrt(sigma2)
  d       <- z^2
  aux     <- signif(lambda*z,digits=7)
  aux1    <- signif(pmax(aux,-37),digits=7)
  cte     <- 2/sqrt(2*pi*sigma2)
  vero    <- cte*exp(-0.5*d)*pnorm(aux1)
 }

 if(family=="stn")
 {
  # Theta=[beta,sigma2,lambda,nu]
  if(missing(X)) {X<-matrix(1,n,1)}
  p       <- ncol(X)
  beta    <- signif(theta[1:p],digits=7)
  sigma2  <- signif(theta[p+1],digits=7)
  lambda  <- signif(theta[p+2],digits=7)
  nu      <- optimize(ftn,c(1,30),y=y,X=X,theta=theta)$minimum
  z       <- (y-X%*%beta)/sqrt(sigma2)
  d       <- z^2
  aux     <- signif(lambda*z,digits=7)
  aux1    <- signif(pmax(aux,-37),digits=7)
  cnugama <- signif(2*gamma((nu+1)/2)/(gamma(nu/2)*sqrt(nu*pi)),digits=7)
  vero    <- signif(cnugama/sqrt(sigma2)*(1+d/nu)^(-(nu+1)/2)*pnorm(aux1),digits=7)
 }

 if(family=="ssl")
 {
  #% =Theta=[beta,sigma^2,lambda,nu]
  n       <- length(y)
  if(missing(X)) {X<-matrix(1,n,1)}
  p       <- ncol(X)
  beta    <- theta[1:p]
  sigma2  <- theta[p+1]
  lambda  <- theta[p+2]
  nu      <- optimize(fsl,c(0.1,20),y,X,theta)$minimum
  z       <- (y-X%*%beta)/sqrt(sigma2)
  d       <- z^2
  aux     <- signif(lambda*z,digits=7)
  aux1    <- signif(pmax(aux,-37),digits=7)
  cte     <- nu/((2*pi*sigma2)^0.5)
  fdps    <- cte*gamma(nu+0.5)*pgamma(1,nu+0.5,scale=2/d)/((d/2)^(nu+0.5))
  fdps[which(z==0)]=cte/(nu+0.5)
  vero    <- 2*fdps*pnorm(aux1)
 }

 if(family=="scn")
 {
  #% Theta=[beta,sigma2,lambda]
  if(missing(X)) {X<-matrix(1,n,1)}
  p       <- ncol(X)
  beta    <- theta[1:p]
  mu      <- X%*%beta
  sigma2  <- theta[p+1]
  lambda  <- theta[p+2]
  L1      <- rbind(1e-6,1e-6)
  L2      <- rbind(0.999999,0.999999)
  nugama0 <- matrix(c(0.5,0.5),2,1)
  nugama  <- optim(nugama0,fcn,method='L-BFGS-B',lower=L1,upper=L2,control=list(maxit=400),y=y,X=X,theta=theta)$par
  nu      <- nugama[1]
  gama    <- nugama[2]
  z       <- (y-mu)/sqrt(sigma2)
  aux     <- signif(lambda*z,digits=7)
  aux1    <- signif(pmax(aux,-37),digits=7)
  vero    <- signif(2*(nu*dnorm(y,mu,sqrt(sigma2/gama))+(1-nu)*dnorm(y,mu,sqrt(sigma2)))*pnorm(aux1),digits=7)
 }

 if(family=="sep")
 {
  #% Theta=[beta,sigma2,lambda]
  p       <- ncol(X)
  beta    <- theta[1:p]
  sigma2  <- theta[p+1]
  lambda  <- theta[p+2]
  nu      <- optimize(fep,c(0.55,1),tol=1e-6,y,X,theta)$minimum
  z       <- (y-X%*%beta)/sqrt(sigma2)
  d       <- z^2
  aux     <- signif(lambda*z,digits=7)
  aux1    <- signif(pmax(aux,-37),digits=7)
  cnu     <- 2*nu/(2^(1/(2*nu))*sqrt(sigma2)*gamma(1/(2*nu)))
  vero    <- cnu*exp(-0.5*d^nu)*pnorm(aux1)
 }

 return(-sum(log(vero)))
}


################################################################################
# Simulated envelope

envel <- function(y,X,theta,family="sn",alpha=0.05){
  n=length(y)
  p=ncol(X)
  mu=X%*%theta[1:p]
  d2=as.numeric((y-mu)^2/theta[p+1])
  d2s=sort(d2)
  d2s=t(d2s)
  xq2 <- qchisq(ppoints(n), p)
  replic=200
  Xsim<-matrix(0,replic,n)
  if (family=="sn"){
    for(i in 1:replic) Xsim[i,]<-rchisq(n, 1)
  }
  if (family=="stn") {
    nu=theta[p+3]
    for(i in 1:replic) Xsim[i,]<-rf(n, 1,nu)
  }
  if (family=="ssl"){
    nu=theta[p+3]
    for(i in 1:replic){
      aux=GeraSlash(n,nu)
      Xsim[i,]<-aux^2}
  }
  if (family=="scn"){
    nu=theta[p+3]
    gama=theta[p+4]
    for(i in 1:replic){
      u1<-rchisq(n, 1)/gama
      u2<-rchisq(n, 1)
      u3<-runif(n)
      id<-(u3<nu)
      u2[id]<-u1[id]
      Xsim[i,]=u2
    }
  }
  if (family=="sep"){
    for(i in 1:replic){
      aux=GeraEP(n,nu)
      Xsim[i,]<-aux^2}
  }

  Xsim2<-apply(Xsim,1,sort)
  d21<-rep(0,n)
  d22<-rep(0,n)
  for(i in 1:n){
    d21[i]  <- quantile(Xsim2[i,],alpha/2)
    d22[i]  <- quantile(Xsim2[i,],1-alpha/2)}
  d2med <-apply(Xsim2,1,mean)
  fy <- range(d2s,d21,d22)
  plot(xq2,d2s,xlab = expression(paste("Theoretical ",chi[1]^2, " quantiles")),
       ylab="Sample values and simulated envelope",pch=20,ylim=fy)
  par(new=T)
  plot(xq2,d21,type="l",ylim=fy,xlab="",ylab="")
  par(new=T)
  plot(xq2,d2med,type="l",ylim=fy,xlab="",ylab="",lty="dashed")
  par(new=T)
  plot(xq2,d22,type="l",ylim=fy,xlab="",ylab="")
}
#theta = c(fit.ssmn$beta,fit.ssmn$sigma2,fit.ssmn$lambda)
#envel(y,x1,theta,family="sn",alpha=0.05)


GeraSlash<- function(n,nu){
  u1 <- runif(n)
  u2 <- u1^(1/(nu))
  u3<-rep(0,n)
  for(j in 1:n){u3[j] <-  rnorm(1, 0, sd=1/sqrt(u2[j]))}
  return(u3)
}

GeraEP<- function(n, nu)
{
  t=rgamma(n,0.5*nu,0.5)
  R=t^(1/(2*nu))
  u=sign(rnorm(n))
  y=u*R
  return(y)
}


