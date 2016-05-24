# Birnbaum Saunders Skew Normal density function

dbssn   <- function(ti,alpha=0.5,beta=1,lambda=1.5)
{
 if(!is.numeric(alpha)||!is.numeric(beta)||!is.numeric(lambda))
 if(alpha<=0){stop("alpha must be positive")}
 if(beta<=0) {stop("beta must be positive")}
 at     <- (1/alpha)*(sqrt(ti/beta)-sqrt(beta/ti))
 At     <- (ti^(-1.5)*(ti+beta))/(2*alpha*sqrt(beta))
 pdf    <- 2*dnorm(at)*At*pnorm(lambda*at)
 return(pdf)
}


# Birnbaum Saunders Skew Normal cumulative distribution function
pbssn    <- function(q,alpha=0.5,beta=1,lambda=1.5)
{
  if(!is.numeric(alpha)||!is.numeric(beta)||!is.numeric(lambda))
  if(alpha<=0){ stop("alpha must be positive")}
  if(beta<=0) { stop("beta must be positive") }
  I      <- vector(mode = "numeric", length = length(q))
  for(i in 1:length(q))
  {
    pdf  <- function(x) 2*dnorm((1/alpha)*(sqrt(x/beta)-sqrt(beta/x)))*(x^(-1.5)*(x+beta))/(2*alpha*sqrt(beta))*pnorm(lambda*(1/alpha)*(sqrt(x/beta)-sqrt(beta/x)))
    I[i] <- integrate(pdf,0,q[i])$value
  }
  return(I)
}

qbssn <- function(p,alpha=0.5,beta=1,lambda=1.5)
{
  if(!is.numeric(alpha)||!is.numeric(beta)||!is.numeric(lambda))
  if(alpha<=0){ stop("alpha must be positive")}
  if(beta<=0) { stop("beta must be positive") }
  zp     <- qsn(p,xi=0,omega=1,lambda) #quantile for the skew-normal distribution
  tp     <- (beta/4)*(alpha*zp + sqrt(alpha^2*zp^2 + 4))^2
  return(tp)
}


rbssn<-function(n,alpha=0.5,beta=1,lambda=1.5)
{
  if(!is.numeric(n)||!is.numeric(alpha)||!is.numeric(beta)||!is.numeric(lambda))
  {stop("non-numeric argument to mathematical function")}
  if(alpha<=0){ stop("alpha must be positive")}
  if(beta<=0) { stop("beta must be positive") }
  z      <- rsn(n,0,1,lambda)
  r      <- beta*((alpha*z*0.5)+sqrt((alpha*z*0.5)^2+1))^2
  return(r)
}

#---------------------------------#
# Log-likelihood function of bssn #
#---------------------------------#

logLikbssn <- function(x,alpha=1, beta=1, lambda=0)
{
  return(sum(log(dbssn(x, alpha = alpha,
                       beta     = beta,
                       lambda    = lambda))))
}
#logLikbssn(c(1,2,4,5),alpha=1, beta=3, delta=2)


#--------------------------#
# BS-SN reliability function  #
#--------------------------#

Rebssn <- function(ti,alpha=0.5,beta=1,lambda=1.5)
{
  if(!is.numeric(alpha)||!is.numeric(beta)||!is.numeric(lambda))
  if(alpha<=0){ stop("alpha must be positive")}
  if(beta<=0) { stop("beta must be positive") }
  at     <- (1/alpha)*(sqrt(ti/beta)-sqrt(beta/ti))
  S      <- 1-psn(at,0,1,lambda)
  return(S)
}

#y<-seq(0,2,0.01)
#f1<-Rebssn(y,0.75,1,1)
#f2<-Rebssn(y,1,1,1)
#f3<-Rebssn(y,1.5,1,1)
#f4<-Rebssn(y,2,1,1)
#den<-cbind(f1,f2,f3,f4)

#matplot(y,den,type="l",col=c("deepskyblue4","firebrick1","darkmagenta","aquamarine4"),ylab="S(t)",xlab="t",lwd=2)
#legend(1.5,1,c(expression(alpha==0.75),expression(alpha==1),expression(alpha==1.5),expression(alpha==2)),col=c("deepskyblue4","firebrick1","darkmagenta","aquamarine4"),lty=1:4,lwd=2,seg.len=2,cex=0.9,box.lty=1,bg=NULL)


#--------------------------#
# BS-SN hazard function  #
#--------------------------#


Fbssn <- function(ti,alpha=0.5,beta=1,lambda=1.5)
{
  if(!is.numeric(alpha)||!is.numeric(beta)||!is.numeric(lambda))
  if(alpha<=0){stop("alpha must be positive")}
  if(beta<=0) {stop("beta must be positive") }
  at     <- (1/alpha)*(sqrt(ti/beta)-sqrt(beta/ti))
  At     <- (ti^(-1.5)*(ti+beta))/(2*alpha*sqrt(beta))
  f      <- 2*dnorm(at)*At*pnorm(lambda*at)
  S      <- 1-psn(at,0,1,lambda)
  R      <- f/S
  return(R)
}

#y<-seq(0,2,0.01)
#f1<-Fbssn(y,0.5,1,-1)
#f2<-Fbssn(y,0.5,1,-2)
#f3<-Fbssn(y,0.5,1,-3)
#f4<-Fbssn(y,0.5,1,-4)
#den<-cbind(f1,f2,f3,f4)
#matplot(y,den,type="l",col=c("deepskyblue4","firebrick1","darkmagenta","aquamarine4"),ylab="h(t)",xlab="t",lwd=2)
#legend(0.1,23,c(expression(lambda==0.75),expression(lambda==1),expression(lambda==1.5),expression(lambda==2)),col=c("deepskyblue4","firebrick1","darkmagenta","aquamarine4"),lty=1:4,lwd=2,seg.len=2,cex=0.9,box.lty=1,bg=NULL)


meanbssn <- function(alpha=0.5,beta=1,lambda=1.5)
{
  if(!is.numeric(alpha)||!is.numeric(beta)||!is.numeric(lambda))
  if(alpha<=0){stop("alpha must be positive")}
  if(beta<=0) {stop("beta must be positive") }
  k      <- 1
  W      <- vector(mode = "numeric", length = length(k))
  for(i in 1:length(k))
  {
    f<-function (z) z^k[i]*sqrt(alpha^2*z^2 + 4)*dsn(z,0,1,lambda)
    W[i]<-integrate(f,lower = -Inf, upper = Inf)$value
  }

  meanbssnresult = (beta/2)*(2+alpha^2+alpha*W[1])
  return(meanbssnresult)
}


varbssn    <- function(alpha=0.5,beta=1,lambda=1.5)
{
  if(!is.numeric(alpha)||!is.numeric(beta)||!is.numeric(lambda))
  if(alpha<=0){stop("alpha must be positive")}
  if(beta<=0) {stop("beta must be positive") }
  k        <- seq(1,3,2)
  W        <- vector(mode = "numeric", length = length(k))

  for(i in 1:length(k))
  {
    f      <-function (z) z^k[i]*sqrt(alpha^2*z^2 + 4)*dsn(z,0,1,lambda)
    W[i]   <-integrate(f,lower = -Inf, upper = Inf)$value
  }

  ET       <- (beta/2)*(2+alpha^2+alpha*W[1]) #Expectation
  ET2      <- (beta^2/2)*(2+4*alpha^2+3*alpha^4+2*alpha*W[1]+alpha^3*W[2])

  varbssnresult  <- ET2-(ET)^2 #Variance
  return(varbssnresult)
}

skewbssn   <- function(alpha=0.5,beta=1,lambda=1.5)
{
  if(!is.numeric(alpha)||!is.numeric(beta)||!is.numeric(lambda))
  if(alpha<=0){stop("alpha must be positive")}
  if(beta<=0) {stop("beta must be positive") }
  k        <- seq(1,7,2)
  W        <- vector(mode = "numeric", length = length(k))
  for(i in 1:length(k))
  {
    f      <-function (z) z^k[i]*sqrt(alpha^2*z^2 + 4)*dsn(z,0,1,lambda)
    W[i]   <-integrate(f,lower = -Inf, upper = Inf)$value
  }

  ET       <- (beta/2)*(2+alpha^2+alpha*W[1]) #Expectation
  ET2      <- (beta^2/2)*(2+4*alpha^2+3*alpha^4+2*alpha*W[1]+alpha^3*W[2])
  ET3      <- (beta^3/2)*(2+9*alpha^2+18*alpha^4+15*alpha^6+3*alpha*W[1]+4*alpha^3*W[2]+alpha^5*W[3])
  var      <- ET2-(ET)^2 #Variance
  skewbssnresult <- (ET3-3*ET*ET2+2*(ET)^3)/(var)^1.5 #skewness
  return(skewbssnresult)
}
#skewbssn(alpha=0.5,beta=1,lambda=1.5)

kurtbssn <- function(alpha=0.5,beta=1,lambda=1.5)
{
  if(!is.numeric(alpha)||!is.numeric(beta)||!is.numeric(lambda))
  if(alpha<=0){stop("alpha must be positive")}
  if(beta<=0) {stop("beta must be positive") }
  k        <- seq(1,7,2)
  W        <- vector(mode = "numeric", length = length(k))

  for(i in 1:length(k))
  {
    f      <- function (z) z^k[i]*sqrt(alpha^2*z^2 + 4)*dsn(z,0,1,lambda)
    W[i]   <- integrate(f,lower = -Inf, upper = Inf)$value
  }

  ET       <- (beta/2)*(2+alpha^2+alpha*W[1])
  ET2      <- (beta^2/2)*(2+4*alpha^2+3*alpha^4+2*alpha*W[1]+alpha^3*W[2])
  ET3      <- (beta^3/2)*(2+9*alpha^2+18*alpha^4+15*alpha^6+3*alpha*W[1]+4*alpha^3*W[2]+alpha^5*W[3])
  ET4      <- (beta^4/2)*(2+16*alpha^2+60*alpha^4+120*alpha^6+105*alpha^8+4*alpha*W[1]+10*alpha^3*W[2]+6*alpha^5*W[3]+alpha^7*W[4])
  var      <- ET2-(ET)^2 #Variance
  kurtbssnresult <- (ET4-4*ET*ET3+6*ET^2*ET2-3*ET^4)/(var)^2 #kurtosis
  return(kurtbssnresult)
}


#Observed Information Matrix

Infmatrix <- function(ti,alpha,beta,lambda)
{
  if(!is.numeric(alpha)||!is.numeric(beta)||!is.numeric(lambda))
  if(alpha<=0){stop("alpha must be positive")}
  if(beta<=0) {stop("beta must be positive") }
  n      <- length(ti)
  at     <- (1/alpha)*(sqrt(ti/beta)-sqrt(beta/ti))
  W      <- dnorm(lambda*at)/pnorm(lambda*at)
  W1     <- -W*(lambda*at+W)
  Iaa    <- (n/alpha^2)+((6*n)/alpha^4)-(3/(alpha^4*beta))*sum(ti)-((3*beta)/alpha^4)*sum(1/ti)+((2*lambda)/alpha^3)*sum((sqrt(ti/beta)-sqrt(beta/ti))*W)+(lambda^2/alpha^4)*sum(((ti/beta)+(beta/ti)-2)*W1)
  Iab    <- -(1/(alpha^3*beta^2))*sum(ti)+(1/alpha^3)*sum(1/ti)+(lambda/(2*alpha^2))*sum(((ti^0.5/beta^1.5)+(1/(beta^0.5*ti^0.5)))*W)+(lambda^2/(2*alpha^3))*sum(((ti/beta^2)-(1/ti))*W1)
  Ibb    <- (n/(2*beta^2))-sum(1/(ti+beta)^2)-(1/(alpha^2*beta^3))*sum(ti)+(lambda/(4*alpha))*sum((((3*ti^0.5)/beta^2.5)+(1/(beta^1.5*ti^0.5)))*W)+(lambda^2/(4*alpha^2))*sum(((ti/beta^3)+(2/beta^2)+(1/(beta*ti)))*W1)
  Ial    <- -(1/alpha^2)*sum((sqrt(ti/beta)-sqrt(beta/ti))*W)-(lambda/alpha^3)*sum(((ti/beta)+(beta/ti)-2)*W1)
  Ibl    <- -(0.5/alpha)*sum(((ti^0.5/beta^1.5)+(1/(beta^0.5*ti^0.5)))*W)+(lambda/(2*alpha^2))*sum(((1/ti)-(ti/beta^2))*W1)
  Ill    <- (1/alpha^2)*sum(((ti/beta)+(beta/ti)-2)*W1)
  r1     <- cbind(Iaa,Iab,Ial)
  r2     <- cbind(Iab,Ibb,Ibl)
  r3     <- cbind(Ial,Ibl,Ill)
  Iobs   <- rbind(r1,r2,r3)
  Iobs   <- as.matrix(-Iobs)
  return(Iobs)
}
#Infmatrix(ti,alpha,beta,delta)

#For computer initial values of beta

mmmeth <- function(ti)
{
  S <- function(ti)
  {
    n <- length(ti)
    return((1/n)*sum(ti))
  }

  R <- function(ti)
  {
    n <- length(ti)
    return(((1/n)*sum(ti^(-1)))^(-1))
  }

  beta0ini  <- (S(ti)*R(ti))^0.5
  alpha0ini <- sqrt(2)*((S(ti)/R(ti))^0.5 - 1)^0.5

  result   <- list(beta0init = beta0ini,alpha0ini=alpha0ini, n = length(ti))
  return(result)
}



#Mixture functions

initialmixbs  <-function(mu=c(1,1),sigma=c(2,2))
{
  g      <-length(mu)
  beta   <-matrix(0,nrow=g,ncol=1)
  alpha  <-matrix(0,nrow=g,ncol=1)
  for(i in 1:g)
  {
    b        <- 0.5*(3*mu[i]-5)
    c        <- 0.5*(5*mu[i]-5*mu[i]^2+sigma[i])
    beta[i]  <- (-b + sqrt(b^2 - 4 * c)) / 2
    alpha[i] <- sqrt(2*(mu[i]/beta[i]-1))
  }
  return(cbind(alpha,beta))
}

#initialmixbs(mu=c(5,5),sigma=c(1,1))


#------------------------------------------------#
#  Random numbers generator function for mixture #
#------------------------------------------------#

rmixbssn <- function(n,alpha,beta,lambda,pii)
{
  y <- vector()
  G <- length(alpha)

  z <- sample(G,size=n,replace=TRUE,prob=pii)
  for(i in 1:n)
  {
    y[i] = rbssn(1,alpha[z[i]],beta[z[i]],lambda[z[i]])
  }
  return(list(z=z,y=y))
}

#rmixbssn(50,alpha=c(0.75,0.25),beta=c(1,1.5),lambda=c(3,2),c(0.8,0.2))




#------------------------------------------------#
#  Density function for mixture BSSN             #
#------------------------------------------------#

d.mixed.bssn <- function(ti, pii, alpha, beta, lambda)
{
  # x: vetor de dados
  # outros parametros devem ser do tipo vetor c() de dimensao g (qtd de misturas)
  g    <- length(pii)
  dens <- 0
  for (j in 1:g) dens <- dens + pii[j]*dbssn(ti, alpha[j], beta[j], lambda[j])
    return(dens)
}

#d.mixed.bssn(ti, pii=c(0.2,0.8), alpha=c(1,1), beta=c(1,1), lambda=c(1.5,2.5))
