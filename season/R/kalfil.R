# Forward and backward sweep of the Kalman filter
`kalfil` <-
function(data,f,vartheta,w,tau,lambda,cmean){
# Setting up matrices
  k<-length(f); # Number of frequencies
  kk<-2*(k+1);
  n<-length(data);
  data<-c(0,data); # Add zero to start of data
  F<-rep(c(1,0),k+1);
  G<-matrix(0,kk,kk);
  G[1,1]=1; G[1,2]=lambda; G[2,2]=1; 
  V<-matrix(0,kk,kk);
  V[1,1]=(tau[1]^2)*(lambda^3)/3; V[1,2]=(tau[1]^2)*(lambda^2)/2; 
  V[2,1]=(tau[1]^2)*(lambda^2)/2; V[2,2]=(tau[1]^2)*lambda; # Trend variance
  for (index in 1:k){
     G[(2*index)+1,(2*index)+1]<-2*cos(2*pi/f[index]); 
     G[(2*index)+1,(2*index)+2]<- -1; G[(2*index)+2,(2*index)+1]<- 1; # Seasonal component
     V[(2*index)+1,(2*index)+1]<-(tau[index+1]^2)*w[index]^2; # w = Seasonal variance, lambda in paper
  }
  C_j<-array(0,c(kk,kk,n+1));
  for (index in 1:kk){
     C_j[index,index,1]<-cmean[index]; # Gives a vague prior for alpha_0
  }
  a_j<-matrix(0,kk,n+1); p_j=matrix(0,kk,n+1);
  e_j<-matrix(0,1,n+1); R_j<-array(0,c(kk,kk,n+1));
# Forward sweep of Kalman filter
  p_j[1,1]<-mean(data[2:(n+1)]); # first obs=mean(p_0);
  ind<-t(0:n);
  for (t in 1:n){
     a_j[,t+1]<-G%*%p_j[,t]; # prediction eqn
     R_j[,,t+1]<-(G%*%C_j[,,t]%*%t(G))+V; # prediction eqn
     e_j[,t+1]<-data[t+1]-(t(F)%*%a_j[,t+1]); # residual
     Q_j<-t(F)%*%R_j[,,t+1]%*%F+(vartheta^2); # fitted value variance
     # Kalman filter:
     p_j[,t+1]<-a_j[,t+1]+(R_j[,,t+1]%*%F%*%(qr.solve(Q_j))%*%e_j[,t+1]);
     C_j[,,t+1]<-R_j[,,t+1]-(R_j[,,t+1]%*%F%*%(qr.solve(Q_j))%*%t(F)%*%t(R_j[,,t+1]));
  }
## Backward sweep of Kalman filter
# DLM matrices
   h_j<-matrix(0,kk,n+1); alpha_j<-matrix(0,kk,n+1);
   B_j<-matrix(0,kk,kk); HH_j<-matrix(0,kk,kk);
# Initial alpha sample;
   alpha_j[,n+1]<-mvrnorm(n = 1, mu=p_j[,n+1], Sigma=C_j[,,n+1]);
# Iterate backwards in time;
   for (t in n:1){
      B_j<-C_j[,,t]%*%t(G)%*%qr.solve(R_j[,,t+1]);
      HH_j<-C_j[,,t]-(B_j%*%R_j[,,t+1]%*%t(B_j));
      h_j[,t]=p_j[,t]+(B_j%*%(alpha_j[,t+1]-a_j[,t+1]));
## Sample alpha_j from a multivariate Normal (mvrnorm from MASS library) 
      alpha_j[,t]<-mvrnorm(n = 1, mu=h_j[,t], Sigma=t(HH_j));
   }
# Estimate vartheta and lambdas (labelled 'w')
   se<-matrix(0,n); # squared error
   alphase<-matrix(0,n-1,k);
   for (t in 2:(n+1)){ #<- time = 1 to n;
      se[t-1]<-(data[t]-(t(F)%*%alpha_j[,t]))^2;
      if (t>2){ # <- 2 to n;
         past<-G%*%alpha_j[,t-1];
         for (index in 1:k){ 
            alphase[t-2,index]<-(alpha_j[(2*index)+1,t]-past[(2*index)+1])^2; 
         }
      }
   }
   shape<-(n/2)-1;
   scale<-sum(se)/2;
   vartheta<-sqrt(rinvgamma(1,shape,scale)) # Update vartheta from inverse gamma;
   shape<-((n-1)/2)-1;
   for (index in 1:k){ 
      scale<-sum(alphase[,index])/2;
# Update lambda from inverse gamma (divide by tau);
      w[index]<-sqrt(rinvgamma(1,shape,scale)/(tau[index+1]^2))
   }
# Estimate amplitude and phase using the periodogram
  amp<-vector(mode="numeric",length=k)
  phase<-vector(mode="numeric",length=k)
  for (index in 1:k){ 
     s<-alpha_j[(2*index)+1,1:n]  
     peri<-peri(s,plot=FALSE)
     loc<-sum(as.numeric(rank(abs(peri$c-f[index]))==1)*(1:length(peri$c))) # Get closest frequency
     amp[index]<-peri$amp[loc]
     phase[index]<-peri$phase[loc]
  }
# Get mean values of C for better starting values
  cmean<-vector(length=kk,mode='numeric')
  for (index in 1:kk){
     cmean[index]<-mean(C_j[index,index,]);
  }  
# Returns
   return(list(vartheta=vartheta,w=w,alpha=alpha_j,amp=amp,phase=phase,cmean=cmean))
} # End of function

