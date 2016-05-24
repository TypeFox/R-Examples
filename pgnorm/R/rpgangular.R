rpgangular <-
function(n,p){

# A function implemented by Steve Kalke

# Description: 
# Samples from the univariate angular distribution corresponding to the p-generalized uniform distribution

# Arguments: 
# p- a positiv constant (default: p=2)
# n- the number of random variables to be simulated

if(missing(p)){p<-2}

V<-matrix(rep(0,2*n),nrow=n)

for (i in 1:n){
U<-runif(2)
while(U[1]^p+U[2]^p>1){U<-runif(2)}
V[i,]<-U }

S1<-sample(0:1,n,replace=T)
S2<-sample(c(-1,1),n,replace=T)
return(S2*atan(V[,1]/V[,2])+pi/2+pi*S1)
}
