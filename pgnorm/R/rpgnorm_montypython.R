rpgnorm_montypython <-
function(n,p){

# A function implemented by Steve Kalke

# Description: 
# Samples from the (univariate, central) p-generalized normal distribution using the Monty Python method

# Arguments: 
# p- a positiv constant (default: p=2)
# n- the number of random variables to be simulated

if (missing(p)){p<-2}

if(p<1){ if( p>=1 ){stop("p-value not valid")}
   else{#data(datasetpgnmp1)
B<-pgnorm::datasetpgnmp1
psi<-B[,2]%*%(B[,1]==p)
beta<-B[,3]%*%(B[,1]==p)
}
}
if(p<=13){#data(datasetpgnmp2)
A<-pgnorm::datasetpgnmp2
b<-A[,2]%*%(A[,1]==(ceiling(p*10^2)*10^(-2)) )
a<-(-p*( log( gamma(1/p) ) + log(1/b) + log( (p^(1/p-1)) ) )  )^(1/p)
}
else{
b<-(gamma(1/p)*(p^(1/p-1)))
a<-0
}
s<-a/(b-a)
R<-rep(0,n)

for (i in 1:n){
U_1<-b*runif(1)
U_2<-1/b*runif(1)

if(U_1<=a) {R[i]<-U_1}

else{
if(U_2<=2*dpgnorm(U_1,p)){R[i]<-U_1}

else{

if(U_2>=1/b-s*( 2*dpgnorm(s*(b-U_1),p) - 1/b ) ){R[i]<-s*(b-U_1)}

else{ U<-runif(2)
if(p>=1){while(U[2]>=(U[1]^(-1))*exp(b^p/p-1/p*(b-log(U[1])/(b^(p-1)))^p)){U<-runif(2)}

R[i]<-b-log(U[1])/(b^(p-1))
}
else{while(U[2]*exp(-b^p/p)*(U[1]^(beta/(beta-1)))>exp( -1/p*(( b+1/psi*( U[1]^(1/(1-beta))-1 )   )^p) ) ){U<-runif(2)}
R[i]<-b+1/psi*( U[1]^(1/(1-beta))-1 )}

}
}

}

}
R2<-sample(c(-1,1),n,replace=T)*R
return(R2)
}
