rpgnorm_ziggurat <-
function(n,p,x){

# A function implemented by Steve Kalke

# Description: 
# Samples from the (univariate, central) p-generalized normal distribution using the Ziggurat method

# Arguments: 
# p- a positiv constant (default: p=2)
# n- the number of random variables to be simulated

if (missing(p)){p<-2}

if(missing(x)){x<-zigsetup(p,2^8)}

if(p<1){       #data(datasetpgnzig)
B<-pgnorm::datasetpgnzig
psi<-B[,2]%*%(B[,1]==p)
beta<-B[,3]%*%(B[,1]==p)
  }
R<-rep(0,n)

for (i in 1:n){
while (R[i]==0){
v<-sample(1:(2^8),1)
if (v==2^8){b<-x[v-1]+(1-pgamma(x[v-1]^p/p,1/p))/2/dpgnorm(x[v-1],p) #x[v-1]*(1-dpgnorm(x[v-2],p)/dpgnorm(x[v-1],p))
u_1<-b*runif(1)
if(u_1<=x[v-1]){R[i]<-u_1}
else{

U<-runif(2)
if(p>=1){
while(U[2]>=(U[1]^(-1))*exp(x[v-1]^p/p-1/p*(x[v-1]-log(U[1])/(x[v-1]^(p-1)))^p)){U<-runif(2)}
R[i]<-x[v-1]-log(U[1])/(x[v-1]^(p-1))
   }
else{
while( U[2]*exp(-x[v-1]^p/p)*(U[1]^(beta/(beta-1)))>exp( -1/p*( ( x[v-1]+1/psi*( U[1]^(1/(1-beta))-1 ))^p) )  ){U<-runif(2)}
R[i]<-x[v-1]+1/psi*( U[1]^(1/(1-beta))-1 )
     }

}

}
else{
u_1<-x[v]*runif(1)
if(v>1){if (u_1<=x[v-1]){R[i]<-u_1}
else{u_2<-2*(dpgnorm(x[v-1],p)-dpgnorm(x[v],p))*runif(1)+2*dpgnorm(x[v],p)

if (u_2<=2*dpgnorm(u_1,p)){R[i]<-u_1}

}
   }
else{ u_2<-2*(dpgnorm(0,p)-dpgnorm(x[1],p))*runif(1)+2*dpgnorm(x[1],p)
if(u_2<=2*dpgnorm(u_1,p)){R[i]<-u_1}
     }


}
}
}

R2<-sample(c(-1,1),n,replace=T)*R
return(R2)

}
