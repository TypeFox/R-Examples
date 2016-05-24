rpgnorm_pgenpolar <-
function(n,p){

# A function implemented by Steve Kalke

# Description: 
# Samples from the (univariate, central) p-generalized normal distribution using the p-generalized polar method

# Arguments: 
# p- a positiv constant (default: p=2)
# n- the number of random variables to be simulated

if (missing(p)){p<-2}

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
   }

V<-rpgunif(n,p)
Y<-V*R
return(Y[,1])
}
