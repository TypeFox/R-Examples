`m4plEstimateMore` <-
function(x,s=1/1.702,b=0,c=0,d=1,m=0,model="T",prior="uniform"){
 w <- 20; NALT <- 5000; p <- 1/NALT; alpha0 <- round(w*p + 1); beta0 <- round(w*(1-p) + 1)
 priori <- function(theta,m,prior) {
  if (prior == "none")    return(1)
  if (prior == "uniform") return(dunif(theta,m-4, m+4))
  if (prior == "normal")  return(dnorm(theta,m,1))
  }
 L    <- function(x,m=0,theta,S=0,C=0,D=0,s=1/1.702,b=0,c=0,d=1,prior) {
                  -log(priori(theta,m,prior) *
                       prod((pm4pl(theta,S=S,C=C,D=D,s=s,b=b,c=c,d=d))^x *
                       (1-pm4pl(theta,S=S,C=C,D=D,s=s,b=b,c=c,d=d))^(1-x)))
  }
 T    <- function(P) L(x=x,theta=P[1],                      m=m,s=s,b=b,c=c,d=d,prior=prior)
 TS   <- function(P) L(x=x,theta=P[1],S=P[2],               m=m,s=s,b=b,c=c,d=d,prior=prior)
 TA   <- function(P) L(x=x,theta=P[1],S=1/P[2],             m=m,s=s,b=b,c=c,d=d,prior=prior)
 TC   <- function(P) L(x=x,theta=P[1],C=P[2],               m=m,s=s,b=b,c=c,d=d,prior=prior)
 TD   <- function(P) L(x=x,theta=P[1],D=P[2],               m=m,s=s,b=b,c=c,d=d,prior=prior)
 TSC  <- function(P) L(x=x,theta=P[1],S=P[2],C=P[3],        m=m,s=s,b=b,c=c,d=d,prior=prior)
 TSD  <- function(P) L(x=x,theta=P[1],S=P[2],D=P[3],        m=m,s=s,b=b,c=c,d=d,prior=prior)
 TCD  <- function(P) L(x=x,theta=P[1],C=P[2],D=P[3],        m=m,s=s,b=b,c=c,d=d,prior=prior)
 TSCD <- function(P) L(x=x,theta=P[1],S=P[2],C=P[3],D=P[4], m=m,s=s,b=b,c=c,d=d,prior=prior)
 res  <- switch(model,
  T    = nlminb(start=0,             objective=T,    lower=-4,           upper=4)$par,
  S    = nlminb(start=c(0,1),        objective=TS,   lower=c(-4,0),      upper=c(4,4))$par,
  A    = nlminb(start=c(0,1),        objective=TA,   lower=c(-4,0),      upper=c(4,4))$par,
  C    = nlminb(start=c(0,0),        objective=TC,   lower=c(-4,0),      upper=c(4,1))$par,
  D    = nlminb(start=c(0,0),        objective=TD,   lower=c(-4,0),      upper=c(4,1))$par,
  SC   = nlminb(start=c(0,1,0),      objective=TSC,  lower=c(-4,0,0),    upper=c(4,4,1))$par,
  SD   = nlminb(start=c(0,1,0),      objective=TSD,  lower=c(-4,0,0),    upper=c(4,4,1))$par,
  CD   = nlminb(start=c(0,0,0),      objective=TCD,  lower=c(-4,0,0),    upper=c(4,1,1))$par,
  SCD  = nlminb(start=c(0,1,0,0),    objective=TSCD, lower=c(-4,0,0,0),  upper=c(4,4,1,1))$par
  )
 if (model != "T") {
  der2       <- fsecond(res, h=.01, FUN=paste("T",model,sep=""))
  adequation <- all(checkAdequation(der2) == TRUE)
  # print(der2); print(checkAdequation(der2))
  if (adequation == TRUE) {
   der2Solve <- solve(der2)
   tse       <- sqrt(diag(der2Solve))
   tcorr     <- cov2cor(der2Solve)
   }
  if (adequation == FALSE) {
   sizeD <- length(res)
   tse   <- rep(NA,sizeD)
   tcorr <- matrix(NA, ncol=sizeD, nrow=sizeD)
   }
  }
 if (model == "T") {
  der2  <- fsecond(res, h=.01, FUN="T")
  if (der2 > 0) tse <- sqrt(1/der2) else tse <- NA
  tcorr <- 1
  }
 LL    <- switch(model,
   T    = -L(x=x,m=m,theta=res,   S=0,     C=0,      D=0,     s=s,b=b,c=c,d=d,prior=prior),
   A    = -L(x=x,m=m,theta=res[1],S=res[2],C=0,      D=0,     s=s,b=b,c=c,d=d,prior=prior),
   S    = -L(x=x,m=m,theta=res[1],S=res[2],C=0,      D=0,     s=s,b=b,c=c,d=d,prior=prior),
   C    = -L(x=x,m=m,theta=res[1],S=0,     C=res[2], D=0,     s=s,b=b,c=c,d=d,prior=prior),
   D    = -L(x=x,m=m,theta=res[1],S=0,     C=0,      D=res[2],s=s,b=b,c=c,d=d,prior=prior),
   SC   = -L(x=x,m=m,theta=res[1],S=res[2],C=res[3], D=0,     s=s,b=b,c=c,d=d,prior=prior),
   SD   = -L(x=x,m=m,theta=res[1],S=res[2],C=0,      D=res[3],s=s,b=b,c=c,d=d,prior=prior),
   CD   = -L(x=x,m=m,theta=res[1],S=0,     C=res[2], D=res[3],s=s,b=b,c=c,d=d,prior=prior),
   SCD  = -L(x=x,m=m,theta=res[1],S=res[2],C=res[3], D=res[4],s=s,b=b,c=c,d=d,prior=prior)
  )

 names(res)  <- switch(model,
  T    = c("T"),
  S    = c("T","S"),
  A    = c("T","A"),
  C    = c("T","C"),
  D    = c("T","D"),
  SC   = c("T","S","C"),
  SD   = c("T","S","D"),
  CD   = c("T","C","D"),
  SCD  = c("T","S","C","D")
  )
 names(LL)  <- model
 names(tse) <- names(res)
 if (length(res) >  1) rownames(tcorr) <- colnames(tcorr) <- names(res)
 if (length(res) == 1) names(tcorr)    <- names(res)
 nParameters   <- length(tse)
 AIC           <-  2*LL - 2*nParameters; BIC <- 2*LL - nParameters*log(length(x))
 logLikelihood <- data.frame(LL,AIC,BIC)
 res           <- list(parameters=res, se=tse, corr=tcorr, logLikelihood=logLikelihood, observeInfo=der2)
 return(res)
 }
 