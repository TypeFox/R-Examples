hmhmm=function(M=100,y){
# Estimation of a hidden Markov model with 2 hidden x 4 observed states

T=length(y)

# Likehood function using the forward/backward algorithm
likej=function(vec){

  P=matrix(c(vec[1],1-vec[1],1-vec[2],vec[2]),ncol=2,byrow=TRUE)
  Q1=vec[3:6]
  Q2=vec[7:10]

  pxy=vec[1:2]/sum(vec[1:2]) # stationary of P

  pyy=rep(1,T)
  pyy[1]=pxy[1]*Q1[y[1]]+pxy[2]*Q2[y[1]]

  for (t in 2:T){

    pxy=pxy[1]*Q1[y[t-1]]*P[1,]+pxy[2]*Q2[y[t-1]]*P[2,]
    pxy=pxy/sum(pxy)
    pyy[t]=(pxy[1]*Q1[y[t]]+pxy[2]*Q2[y[t]])
    }

  sum(log(pyy))
}

propal=function(vec1,vec2){

# Warning : this is the inverse proposal ratio for a vec1 -> vec2 move
# ie log( q(vec2|vec1) / q(vec1|vec2) )
# This is a ratio hence the 8 lgamma(alf) vanish.

   -sum(lgamma(1+(alf*(1-vec1[1:2]))))+sum(lgamma(1+(alf*(1-vec2[1:2]))))-
      sum(lgamma(1+(alf*vec1)))+sum(lgamma(1+(alf*vec2)))+
      sum(log(vec2)*(alf*vec1))-sum(log(vec1)*(alf*vec2))+
      sum(log(1-vec2[1:2])*(alf*(1-vec1[1:2])))-sum(log(1-vec1[1:2])*(alf*(1-vec2[1:2])))
}

# MCMC HM on P and Q1, Q1
# which means 10 parameters (with redundancy)

# First, get a starting point
BigR=matrix(0,ncol=10,nrow=M) #MCMC sample
olike=rep(0.123456789,M)      #corresponding likehood

# Parameters of P
a=.5*runif(1) #0.7  #runif(1)
b=.5*runif(1) #0.58 #runif(1)

# Parameters of Q1
p1=runif(4) #c(0.378,0.0244,0.407, 0.189) 
p1=p1/sum(p1)

# Parameters of Q2
p2=runif(4) #c(0.313,0.403,0.020,0.262)  
p2=p2/sum(p2)

BigR[1,]=c(a,b,p1,p2)
curlike=likej(BigR[1,])
olike[1]=curlike

ace=0 # acceptance rate

for (r in 2:M){

  alf=sample(c(1,10,100,1000,10000,100000),1) # scale of the proposal

  # Dirichlet proposals (the plus 1 are the priors)
  
  na = 1/(1+rgamma(1,1+(alf*(1-BigR[r-1,1])),1)/rgamma(1,1+(alf*(BigR[r-1,1])),1))
  nb = 1/(1+rgamma(1,1+(alf*(1-BigR[r-1,2])),1)/rgamma(1,1+(alf*(BigR[r-1,2])),1))

  p1=rgamma(4,1+(alf*BigR[r-1,3:6]))
  p2=rgamma(4,1+(alf*BigR[r-1,7:10]))
  p1=p1/sum(p1)
  p2=p2/sum(p2)
  
  propa=c(na,nb,p1,p2)

  # Acceptance prob

  proplike=likej(propa)
  accp=proplike-curlike+propal(propa,BigR[r-1,])

  # Decision

  if (log(runif(1))<accp){

     curlike=proplike
     BigR[r,]=propa
     ace = ace+1

     }else{

	BigR[r,]=BigR[r-1,]
     }
  olike[r]=curlike
}

list(par=BigR,olike=olike)
}
