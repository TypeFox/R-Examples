pearson0findM <- function(mean,variance,skewness,kurtosis,moments) {
  if (!missing(moments)) { mmm <- moments[[1]]; vvv <- moments[[2]]; 
                           sss <- moments[[3]]; kkk <- moments[[4]]
  } else {
    mmm <- mean; vvv <- variance; sss <- skewness; kkk <- kurtosis  
  }
  pearsonFitM(mmm,vvv,0,3)
}

pearsonIfindM <- function(mean,variance,skewness,kurtosis,moments) {
  if (!missing(moments)) { mmm <- moments[[1]]; vvv <- moments[[2]]; 
                           sss <- moments[[3]]; kkk <- moments[[4]]
  } else {
    mmm <- mean; vvv <- variance; sss <- skewness; kkk <- kurtosis  
  }
  if (isTRUE(all.equal(kkk,1))||(kkk<1)) {
    kkk <- 2
  }
  slim1 <- max(0,-2+2/3*kkk)
  slim2 <- kkk-1
  if (isTRUE(all.equal(sss^2,slim1))||isTRUE(all.equal(sss^2,slim2))||
      (sss^2>slim2)||(sss^2<slim1)) {
    sss <- sign(sss)*0.5*(sqrt(slim1)+sqrt(slim2))
  }
  pearsonFitM(mmm,vvv,sss,kkk)
}

pearsonIIfindM <- function(mean,variance,skewness,kurtosis,moments) {
  if (!missing(moments)) { mmm <- moments[[1]]; vvv <- moments[[2]]; 
                           sss <- moments[[3]]; kkk <- moments[[4]]
  } else {
    mmm <- mean; vvv <- variance; sss <- skewness; kkk <- kurtosis  
  }
  if (isTRUE(all.equal(kkk,1))||isTRUE(all.equal(kkk,3))||(kkk<1)||(kkk>3)) {
    kkk <- 2
  }
  pearsonFitM(mmm,vvv,0,kkk)
}

pearsonIIIfindM <- function(mean,variance,skewness,kurtosis,moments) {
  if (!missing(moments)) { mmm <- moments[[1]]; vvv <- moments[[2]]; 
                           sss <- moments[[3]]; kkk <- moments[[4]]
  } else {
    mmm <- mean; vvv <- variance; sss <- skewness; kkk <- kurtosis  
  }
  if (isTRUE(all.equal(kkk,3))||(kkk<3)) {
    kkk <- 4
  }
  pearsonFitM(mmm,vvv,sign(sss)*sqrt(-2+2/3*kkk),kkk)
}

pearsonIVfindM <- function(mean,variance,skewness,kurtosis,moments) {
  if (!missing(moments)) { mmm <- moments[[1]]; vvv <- moments[[2]]; 
                           sss <- moments[[3]]; kkk <- moments[[4]]
  } else {
    mmm <- mean; vvv <- variance; sss <- skewness; kkk <- kurtosis  
  }
  if (isTRUE(all.equal(kkk,3))||(kkk<3)) {
    kkk <- 4
  }
  slim2 <- (kkk^2+78*kkk-63-sqrt(kkk+147)*sqrt(kkk+3)^3)/72
  if (isTRUE(all.equal(sss,0))||isTRUE(all.equal(sss^2,slim2))||(sss^2>slim2)) {
    sss <- sign(sss)*0.5*sqrt(slim2)
  }
  pearsonFitM(mmm,vvv,sss,kkk)
}

pearsonVfindM <- function(mean,variance,skewness,kurtosis,moments) {
  if (!missing(moments)) { mmm <- moments[[1]]; vvv <- moments[[2]]; 
                           sss <- moments[[3]]; kkk <- moments[[4]]
  } else {
    mmm <- mean; vvv <- variance; sss <- skewness; kkk <- kurtosis  
  }
  if (isTRUE(all.equal(kkk,3))||(kkk<3)) {
    kkk <- 4
  }
  pearsonFitM(mmm,vvv,
              sign(sss)*sqrt((kkk^2+78*kkk-63-sqrt(kkk+147)*sqrt(kkk+3)^3)/72),
              kkk)
}

pearsonVIfindM <- function(mean,variance,skewness,kurtosis,moments) {
  if (!missing(moments)) { mmm <- moments[[1]]; vvv <- moments[[2]]; 
                           sss <- moments[[3]]; kkk <- moments[[4]]
  } else {
    mmm <- mean; vvv <- variance; sss <- skewness; kkk <- kurtosis  
  }
  if (isTRUE(all.equal(kkk,3))||(kkk<3)) {
    kkk <- 4
  }
  slim1 <- (kkk^2+78*kkk-63-sqrt(kkk+147)*sqrt(kkk+3)^3)/72
  slim2 <- -2+2/3*kkk
  if (isTRUE(all.equal(sss^2,slim1))||isTRUE(all.equal(sss^2,slim2))||
      (sss^2>slim2)||(sss^2<slim1)) {
    sss <- sign(sss)*0.5*(sqrt(slim1)+sqrt(slim2))
  }
  pearsonFitM(mmm,vvv,sss,kkk)
}

pearsonVIIfindM <- function(mean,variance,skewness,kurtosis,moments) {
  if (!missing(moments)) { mmm <- moments[[1]]; vvv <- moments[[2]]; 
                           sss <- moments[[3]]; kkk <- moments[[4]]
  } else {
    mmm <- mean; vvv <- variance; sss <- skewness; kkk <- kurtosis  
  }
  if (isTRUE(all.equal(kkk,3))||(kkk<3)) {
    kkk <- 4
  }
  pearsonFitM(mmm,vvv,0,kkk)
}

empMoments <- function(x) {
  n <- length(x)
  mmm <- mean(x)
  vvv <- var(x)*(n-1)/n
  if (vvv>0) sss <- sum((x-mmm)^3/sqrt(vvv)^3)/n else sss <- 0
  if (vvv>0) kkk <- sum((x-mmm)^4/vvv^2)/n else kkk <- 0
  c(mean=mmm,variance=vvv,skewness=sss,kurtosis=kkk)
}

pearson0fitML <- function(x,...) {                                              # exact MLE available, "fake" nlminb output...
#  sval  <- pearson0findM(moments=empMoments(x))[-1]
#  tfunc <- function(sval) -sum(dpearson0(x,params=sval,log=TRUE))
#  nlminb(sval,tfunc,lower=c(-Inf,0),upper=c(Inf,Inf),...)
  n    <- length(x)
  Mean <- mean(x)
  Sd   <- sd(x)*sqrt((n-1)/n)
  list(par         = c(mean=Mean,sd=Sd),
       objective   = -sum(dnorm(x,mean=Mean,sd=Sd,log=TRUE)),
       convergence = 0, iterations = 1, 
       evaluations = c(`function`=1,gradient=0),
       message     = NULL)
}

pearsonIfitML <- function(x,...) {
  sval  <- pearsonIfindM(moments=empMoments(x))[-1]
  if (sval[[4]]>0) {
    sval[[3]] <- min(sval[[3]],x)-0.1
    sval[[4]] <- max(sval[[4]],x-sval[[3]])+0.1
  } else {
    sval[[3]] <- max(sval[[3]],x)+0.1
    sval[[4]] <- min(sval[[4]],x-sval[[3]])-0.1
  }
  tfunc <- function(sval) {
    tmp <- -sum(dpearsonI(x,params=sval,log=TRUE))
    if (is.na(tmp)) +Inf else tmp
  }  
  nlminb(sval,tfunc,lower=c(0,0,-Inf,-Inf),upper=c(Inf,Inf,Inf,Inf),...)
}

pearsonIIfitML <- function(x,...) {
  sval  <- pearsonIIfindM(moments=empMoments(x))[-1]
  if (sval[[3]]>0) {
    sval[[2]] <- min(sval[[2]],x)-0.01
    sval[[3]] <- max(sval[[3]],x-sval[[2]])+0.01
  } else {
    sval[[2]] <- max(sval[[2]],x)+0.01
    sval[[3]] <- min(sval[[3]],x-sval[[2]])-0.01
  }
  tfunc <- function(sval) {
    tmp <- -sum(dpearsonII(x,params=sval,log=TRUE))
    if (is.na(tmp)) +Inf else tmp
  }      
  nlminb(sval,tfunc,lower=c(0,-Inf,-Inf),upper=c(Inf,Inf,Inf),...)
}

pearsonIIIfitML <- function(x,...) {
  sval  <- pearsonIIIfindM(moments=empMoments(x))[-1]
  if (sval[[3]]>0) sval[[2]] <- min(sval[[2]],x)-0.01 else 
                   sval[[2]] <- max(sval[[2]],x)+0.01
  tfunc <- function(sval) {
    tmp <- -sum(dpearsonIII(x,params=sval,log=TRUE))
    if (is.na(tmp)) +Inf else tmp
  }      
  nlminb(sval,tfunc,lower=c(0,-Inf,-Inf),upper=c(Inf,Inf,Inf),...)
}

pearsonIVfitML <- function(x,...) {
  sval  <- pearsonIVfindM(moments=empMoments(x))[-1]
  tfunc <- function(sval) {
    tmp <- -sum(dpearsonIV(x,params=sval,log=TRUE))
    if (is.na(tmp)) +Inf else tmp
  }      
  nlminb(sval,tfunc,lower=c(0.5,-Inf,-Inf,0),upper=c(Inf,Inf,Inf,Inf),...) 
}

pearsonVfitML <- function(x,...) {
  sval  <- pearsonVfindM(moments=empMoments(x))[-1]
  if (sval[[3]]>0) sval[[2]] <- min(sval[[2]],x)-0.01 else 
                   sval[[2]] <- max(sval[[2]],x)+0.01
  tfunc <- function(sval) { 
    tmp <- -sum(dpearsonV(x,params=sval,log=TRUE))
    if (is.na(tmp)) +Inf else tmp
  } 
  nlminb(sval,tfunc,lower=c(0,-Inf,-Inf),upper=c(Inf,Inf,Inf),...)
}

pearsonVIfitML <- function(x,...) {
  sval  <- pearsonVIfindM(moments=empMoments(x))[-1]
  if (sval[[4]]>0) sval[[3]] <- min(sval[[3]],x)-0.01 else 
                   sval[[3]] <- max(sval[[3]],x)+0.01
  tfunc <- function(sval) {
    tmp <- -sum(dpearsonVI(x,params=sval,log=TRUE))
    if (is.na(tmp)) +Inf else tmp
  } 
  nlminb(sval,tfunc,lower=c(0,0,-Inf,-Inf),upper=c(Inf,Inf,Inf,Inf),...)
}

pearsonVIIfitML <- function(x,...) {
  sval  <- pearsonVIIfindM(moments=empMoments(x))[-1]
  tfunc <- function(sval) {
    tmp <- -sum(dpearsonVII(x,params=sval,log=TRUE))
    if (is.na(tmp)) +Inf else tmp
  } 
  nlminb(sval,tfunc,lower=c(0,-Inf,-Inf),upper=c(Inf,Inf,Inf),...)
}

pearsonFitML <- function(x,...) {
  suffix <- c("0","I","II","III","IV","V","VI","VII")
  res <- vector("list",length(suffix))
  obj <- numeric(length(suffix))
  for (i in seq_along(suffix)) {
    fname  <- paste("pearson",suffix[i],"fitML",sep="")
    res[[i]] <- do.call(fname,args=list(x=x,`...`=...))
    obj[i]   <- res[[i]]$objective
  }
  ML <- which.min(obj)
  as.list(c(type=ML-1,res[[ML]]$par))
}

pearsonFitM <- function(mean,variance,skewness,kurtosis,moments) {
  if (!missing(moments)) { mmm <- moments[[1]]; vvv <- moments[[2]]; 
                           sss <- moments[[3]]; kkk <- moments[[4]]
  } else {
    mmm <- mean; vvv <- variance; sss <- skewness; kkk <- kurtosis  
  }
  if (isTRUE(all.equal(sss^2,kkk-1))) {
    p <- (1+sss/sqrt(4+sss^2))/2
    a <- mmm-vvv*sqrt((1-p)/p)
    b <- mmm+vvv*sqrt(p/(1-p))
    errstr <- paste("Target distribution not in Pearson system,\n",
                  "  try a discrete distribution with mass\n",
                  "  p   = ",p  ," on a = ",a," and mass\n",
                  "  1-p = ",1-p," on b = ",b," instead!",sep="")
    stop(errstr)               
  }
  if (sss^2>kkk-1)
    stop("There are no probability distributions with these moments")
  c0  <- (4*kkk-3*sss^2)/(10*kkk-12*sss^2-18)*vvv
  c1  <- sss*(kkk+3)/(10*kkk-12*sss^2-18)*sqrt(vvv)
  c2  <- (2*kkk-3*sss^2-6)/(10*kkk-12*sss^2-18)
  if (isTRUE(all.equal(sss,0))) {
    if (isTRUE(all.equal(kkk,3))) return(list(type=0,mean=mmm,sd=sqrt(vvv)))    # type 0 (normal distribution)
    if (kkk<3) {                                                                # type II
      a1 <- sqrt(vvv)/2*(-sqrt(-16*kkk*(2*kkk-6))/(2*kkk-6))
      m1 <-  -1/(2*c2)
      stopifnot(m1>-1)
      sca <- 2*a1
      loc <- mmm - sca / 2
      return(list(type=2,a=1+m1,location=loc,scale=sca))                        # -> scale*beta(a,a)+location
    }  
    if (kkk>3) {                                                                # type VII
      r      <- 6*(kkk-1)/(2*kkk-6)
      a      <- sqrt(vvv*(r-1))
      dof    <- 1+r
      return(list(type=7,df=dof,location=mmm,scale=a/sqrt(dof)))                # -> scale*t(df)+location
    }  
  } else if (!isTRUE(all.equal(2*kkk-3*sss^2-6,0))) {
    kap <- 0.25*sss^2*(kkk+3)^2/((4*kkk-3*sss^2)*(2*kkk-3*sss^2-6))
    if (kap<0) {                                                                # type I
      a1 <- sqrt(vvv)/2*((-sss*(kkk+3)-sqrt(sss^2*(kkk+3)^2-4*(4*kkk-3*sss^2)*
            (2*kkk-3*sss^2-6)))/(2*kkk-3*sss^2-6))
      a2 <- sqrt(vvv)/2*((-sss*(kkk+3)+sqrt(sss^2*(kkk+3)^2-4*(4*kkk-3*sss^2)*
            (2*kkk-3*sss^2-6)))/(2*kkk-3*sss^2-6))
      if (a1>0) { tmp <- a1; a1 <- a2; a2 <- tmp }
      a  <- c1
      m1 <- -(sss*(kkk+3)+a1*(10*kkk-12*sss^2-18)/sqrt(vvv))/
             (sqrt(sss^2*(kkk+3)^2-4*(4*kkk-3*sss^2)*(2*kkk-3*sss^2-6)))
      m2 <- -(-sss*(kkk+3)-a2*(10*kkk-12*sss^2-18)/sqrt(vvv))/
             (sqrt(sss^2*(kkk+3)^2-4*(4*kkk-3*sss^2)*(2*kkk-3*sss^2-6)))
      stopifnot((m1>-1)&&(m2>-1))
      sca <- (a2-a1)
      loc <- mmm - sca * (m1+1)/(m1+m2+2)
      return(list(type=1,a=1+m1,b=1+m2,location=loc,scale=sca))                 # -> scale*beta(a,b)+location
    }  
    if (isTRUE(all.equal(kap,1))) {                                             # type V
      C1 <- c1/(2*c2)
      sca <- -(c1-C1)/c2
      return(list(type=5,shape=1/c2-1,location=mmm-C1,scale=sca))               # -> invgamma(shape,scale)+location, TAKE CARE: scale maybe negativ!!!
    }  
    if (kap>1) {                                                                # type VI
      a1 <- sqrt(vvv)/2*((-sss*(kkk+3)-sqrt(sss^2*(kkk+3)^2-4*(4*kkk-3*sss^2)*
            (2*kkk-3*sss^2-6)))/(2*kkk-3*sss^2-6))
      a2 <- sqrt(vvv)/2*((-sss*(kkk+3)+sqrt(sss^2*(kkk+3)^2-4*(4*kkk-3*sss^2)*
            (2*kkk-3*sss^2-6)))/(2*kkk-3*sss^2-6))
      if (a1>0) { tmp <- a1; a1 <- a2; a2 <- tmp }
      a  <- c1
      m1 <- (a + a1)/(c2*(a2-a1))
      m2 <- - (a + a2)/(c2*(a2-a1))
      sca <- (a2-a1)
      loc <- mmm + sca * (m2+1)/(m2+m1+2)
      return(list(type=6,a=1+m2,b=-m2-m1-1,location=loc,scale=sca))             # -> scale*betaprime(a,b)+location
    }  
    if ((kap>0)&&(kap<1)) {                                                     # type IV
      r        <- 6*(kkk-sss^2-1)/(2*kkk-3*sss^2-6)
      nu       <- -r*(r-2)*sss/sqrt(16*(r-1)-sss^2*(r-2)^2)
      scale    <- sqrt(vvv*(16*(r-1)-sss^2*(r-2)^2))/4
      location <- mmm - ((r-2)*sss*sqrt(vvv))/4
      m      <- 1+r/2
      return(list(type=4,m=m,nu=nu,location=location,scale=scale))
    }
  } else {                                                                      # type III
    m <- c0/(c1^2) - 1
    loc <-  mmm - c0/c1
    sca <-  c1
    return(list(type=3,shape=m+1,location=loc,scale=c1))                        # -> gamma(shape,scale)+location, TAKE CARE: scale maybe negativ!!!
  } 
}

pearsonMSC <- function(x,...) {
  suffix <- c("0","I","II","III","IV","V","VI","VII")                           # names of distribution sub-classes
  K      <- c(2,4,3,3,4,3,4,3)                                                  # number of estimated parameters
  names(K) <- suffix
  n      <- length(x)                                                           # number of observations
  res    <- vector("list",length(suffix))
  lnL    <- numeric(length(suffix))
  names(lnL) <- suffix
  for (i in seq_along(suffix)) {
    fname    <- paste("pearson",suffix[i],"fitML",sep="")
    res[[i]] <- do.call(fname,args=list(x=x,`...`=...))
    lnL[i]   <- -res[[i]]$objective
  }
  AIC  <- 2*K-2*lnL
  BIC  <- K*log(n)-2*lnL
  AICc <- 2*K*n/(n-K-1)-2*lnL
  HQ   <- 2*K*log(log(n))-2*lnL
  fits <- lapply(res,function(x) x$par)
  names(fits) <- suffix
  for (i in seq_along(fits)) fits[[i]] <- c(type=i-1,fits[[i]])
  MSCs <- rbind("ML"=-2*lnL,AIC=AIC,AICc=AICc,BIC=BIC,HQC=HQ)
  Best <- vector("list",nrow(MSCs))
  names(Best) <- rownames(MSCs)
  for (i in seq_along(Best)) Best[[i]] <- fits[[which.min(MSCs[i,])]]
  list(MSCs=MSCs,logLik=lnL,FittedDistributions=fits,Best=Best)
}
