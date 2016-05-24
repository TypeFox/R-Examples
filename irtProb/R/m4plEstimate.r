`m4plEstimate` <-
function(x,s=1/1.702,b=0,c=0,d=1,m=0,model="T",prior="uniform"){
 w <- 20; NALT <- 5000; p <- 1/NALT; alpha0 <- round(w*p + 1) ;beta0 <- round(w*(1-p) + 1)
 priori <- function(theta,m,prior) {
  if (prior == "none")    return(1)
  if (prior == "uniform") return(dunif(theta,m-4, m+4))
  if (prior == "normal")  return(dnorm(theta,m,1))
  }
 L    <- function(x,m=0,theta,S=0,C=0,D=0,s=1/1.702,b=0,c=0,d=1,prior) {
                  -log(priori(theta,m,prior) *
                       prod((pm4pl(theta,S=S,C=C,D=D,s=s,b=b,c=c,d=d))^x *
                       (1-pm4pl(theta,S=S,C=C,D=D,s=s,b=b,c=c,d=d))^(1-x))
                       )
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
 return(res)
 }

