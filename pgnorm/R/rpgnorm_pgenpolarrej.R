rpgnorm_pgenpolarrej <-
function(n,p){

# A function implemented by Steve Kalke

# Description: 
# Samples from the (univariate, central) p-generalized normal distribution using the p-generalized rejecting polar method

# Arguments: 
# p- a positiv constant (default: p=2)
# n- the number of random variables to be simulated

if (missing(p)){p<-2}

V<-matrix(rep(0,2*n),nrow=n)

k<-floor(2/p)
p1<-2/p-k
p2<-1-p1
R<-rep(0,n)
for (i in 1:n){
s<-sum(log(runif(k)))
if(k!=2/p){
U<-runif(2)
while (U[1]^(1/p1)+U[2]^(1/p2)>1){U<-runif(2)}
s<-s+log(runif(1))*(U[1]^(1/p1))/ ( U[1]^(1/p1) + U[2]^(1/p2)  )
     }
R[i]<-(-p*s)^(1/p)
U<-runif(2)
while(U[1]^p+U[2]^p>1){U<-runif(2)}
V[i,]<-U
   }

S1<-sample(c(-1,1),n,replace=T)
S2<-sample(c(-1,1),n,replace=T)

U<-cbind(S1*V[,1],S2*V[,2])/((abs(V[,1])^p+abs(V[,2])^p)^(1/p))
Y<-U*R
return(Y[,1])

}
