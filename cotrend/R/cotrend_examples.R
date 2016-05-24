######################################################################
#       		cotremd_examples.R
#
# Examples from paper by Guo and Shintani, version: Feb 2011
# All examples are simple Monte Carlo simulations where performance
# was not major concern.
#
# Examples are numbered following the equation numbers in the paper.
# Not all examples in the paper have been included.
#
# For all examples the length of the time series is given by paramter
# "T". The corect rank is given if with P->1 if T-> inf
#
# Sept,2011
######################################################################

#######################################################################
# pg. 5 equation(3)
# r1=0,r2=1
#######################################################################

example_eq3 <- function(T=50,mu1=0.5,mu2=2,c1=0.5,c2=1){

t        <- seq(1,T,by=1)
epsilon1 <- rnorm(T)
epsilon2 <- rnorm(T)

d1    <- c1+mu1*t
d2    <- c2+mu2*t
s1    <- cumsum(epsilon1)
s2    <- cumsum(epsilon2)
y1    <- d1+s1
y2    <- d2+s2

return(cbind(y1,y2)) # dim(Y)=(T,2)
}

#######################################################################
# pg. 5 equation(4)
# r1=1,r2=1
#######################################################################

example_eq4 <- function(T=50,mu1=0.5,mu2=2,c1=0.5,c2=1){

t        <- seq(1,T,by=1)
epsilon1 <- rnorm(T)
epsilon2 <- rnorm(T)

d1    <- c1+mu1*t
d2    <- c2+mu2*t
s1    <- cumsum(epsilon1)
s2    <- (mu2/mu1)*s1 + epsilon2
y1    <- d1+s1
y2    <- d2+s2

return(cbind(y1,y2)) # dim(Y)=(T,2)
}
#########################################################################
# p. 19 equation (8)
# Default parameters extrated from paper. Rank changes by changing
# rho1.
# rho1 = 1: r1=0, r2 =1 ; rho1 = 0.5 r1=1,r2=1
#########################################################################

example_eq8 <- function(T=50,rho1=1,mu0=2,mu1=0.5,c=0.5,tau=0.5,y0=0){
	
t        <- seq(1,T,by=1)
epsilon1 <- rnorm(T)
yold1    <- y0
Y <- matrix(0,length(t),3)
for(i in 1:length(t)){
  y1    <- rho1*yold1+epsilon1[i]
  yold1 <- y1
  y2    <- c+mu0*t[i]
  if(t[i]<=tau*T){
    y3 <- c+mu0*t[i]  
  }else{
    y3 <- c+(mu0-mu1)*tau*T+mu1*t[i]
  }
  Y[i,]<-c(y1,y2,y3)
}

return(Y) # matrix of dim(Y) = (T,3)

}

#########################################################################
# pg 19: r1= 0 r2 =0, equation(9)
#########################################################################

example_eq9 <- function(T=400,mu0=2,mu1=0.5,c=0.5,tau1=0.5,tau2=1/3){

t        <- seq(1,T,by=1)
epsilon1 <- rnorm(T)
Y <- matrix(0,length(t),3)
for(i in 1:length(t)){
  y1    <- c+mu0*t[i]+epsilon1[i]

  if(t[i]<=tau1*T){
    y2 <- c+mu0*t[i]  
  }else{
    y2 <- c+(mu0-mu1)*tau1*T+mu1*t[i]
  }

  if(t[i]<=tau2*T){
    y3 <- c+mu0*t[i]  
  }else{
    y3 <- c+(mu0-mu1)*tau2*T+mu1*t[i]
  }
  Y[i,] <- c(y1,y2,y3)
}

return(Y) # dim(Y) = (T,3)

}

