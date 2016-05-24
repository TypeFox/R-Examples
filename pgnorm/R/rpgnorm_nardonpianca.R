rpgnorm_nardonpianca <-
function(n,p){

# A function implemented by Steve Kalke

# Description: 
# Samples from the (univariate, central) p-generalized normal distribution using the method of Nardon and Pianca

# Arguments: 
# p- a positiv constant (default: p=2)
# n- the number of random variables to be simulated

if (missing(p)){p<-2}

l<-floor(1/p)
p1<-1/p-l
p2<-1-p1

R<-rep(0,n)

for (i in 1:n){
s<-sum(log(runif(l)))
if(l!=1/p){
U<-runif(2)
while (U[1]^(1/p1)+U[2]^(1/p2)>1){U<-runif(2)}
s<-s+log(runif(1))*(U[1]^(1/p1))/ ( U[1]^(1/p1) + U[2]^(1/p2) )
     }

R[i]<-(-p*s)^(1/p)

}

return(R*sample(c(-1,1),n,replace=TRUE))

}
