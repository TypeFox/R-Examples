# nscosinor.initial.R
# Estimate starting values for seasonal variance based on a linear trend and given values for tau
# Dec 2011

# set up variables and matrices
nscosinor.initial=function(data,response,tau,lambda=1/12,n.season){
k<-1; # Assume just one season
kk<-2*(k+1);
n<-nrow(data);
F<-rep(c(1,0),k+1);
alpha_j<-matrix(0,kk,n+1);
G<-matrix(0,kk,kk);
G[1,1]=1; G[1,2]=lambda; G[2,2]=1; 
# linear model
time=1:n
vector.response=subset(data,select=response)[,1] # replaced attach
model=glm(vector.response~time)
## put predictions into alpha
# 1. trend
alpha_j[1,2:(n+1)]=fitted(model)
# 2. season
sdr=sd(resid(model))
alpha_j[3,2:(n+1)]=(sdr/10)*cos(2*pi*(1:n+1)/12) # sinusoid with amplitude equal to 10% of standard deviation of residuals
# estimate initial value for w
   se<-matrix(0,n); # squared error
   alphase<-matrix(0,n-1,k);
   for (t in 2:(n+1)){ #<- time = 1 to n;
      se[t-1]<-(vector.response[t-1]-(t(F)%*%alpha_j[,t]))^2; 
      if (t>2){ # <- 2 to n;
         past<-G%*%alpha_j[,t-1];
         for (index in 1:k){ 
            alphase[t-2,1]<-(alpha_j[(2*1)+1,t]-past[(2*1)+1])^2; 
         }
      }
   }
# variance initial value
   shape<-(n/2)-1;
   scale<-sum(se)/2;
   vartheta = scale/(shape-1) # based on mean

# seasonal initial values
   shape<-((n-1)/2)-1;
   scale<-sum(alphase)/2;
   w=(scale/(shape-1))/(tau[2]^2) # use mean rather than random sampling
   initial.value=c(vartheta,rep(w,n.season))
   return(initial.value)
}   
