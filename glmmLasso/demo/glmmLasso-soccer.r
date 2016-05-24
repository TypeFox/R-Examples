library(glmmLasso)
data("soccer")
## generalized additive mixed model
## grid for the smoothing parameter

## center all metric variables so that also the starting values with glmmPQL are in the correct scaling

soccer[,c(4,5,9:16)]<-scale(soccer[,c(4,5,9:16)],center=T,scale=T)
soccer<-data.frame(soccer)


lambda <- seq(500,0,by=-5)

family = poisson(link = log)


################## First Simple Method ############################################
## Using BIC (or AIC, respectively) to determine the optimal tuning parameter lambda

BIC_vec<-rep(Inf,length(lambda))

## first fit good starting model
library(MASS);library(nlme)
PQL<-glmmPQL(points~1,random = ~1|team,family=family,data=soccer)
Delta.start<-c(as.numeric(PQL$coef$fixed),rep(0,6),as.numeric(t(PQL$coef$random$team)))
Q.start<-as.numeric(VarCorr(PQL)[1,1])

for(j in 1:length(lambda))
{
print(paste("Iteration ", j,sep=""))
  
glm1 <- try(glmmLasso(points~transfer.spendings  
        + ave.unfair.score + ball.possession
        + tackles + ave.attend + sold.out, rnd = list(team=~1),  
        family = family, data = soccer, lambda=lambda[j],switch.NR=T,final.re=TRUE,
        control=list(start=Delta.start,q_start=Q.start)), silent=TRUE)  


if(class(glm1)!="try-error")
{  
BIC_vec[j]<-glm1$bic
}
        
}
    
opt<-which.min(BIC_vec)
        
glm1_final <- glmmLasso(points~transfer.spendings  
        + ave.unfair.score + ball.possession
        + tackles + ave.attend + sold.out, rnd = list(team=~1),  
        family = family, data = soccer, lambda=lambda[opt],switch.NR=F,final.re=TRUE,
        control=list(start=Delta.start,q_start=Q.start))
         
  
        
summary(glm1_final)





################## Second Simple Method ###########################
## Using 5-fold CV to determine the optimal tuning parameter lambda

### set seed
set.seed(123)
N<-dim(soccer)[1]
ind<-sample(N,N)
lambda <- seq(500,0,by=-5)

kk<-5
nk <- floor(N/kk)

Devianz_ma<-matrix(Inf,ncol=kk,nrow=length(lambda))

## first fit good starting model
library(MASS);library(nlme)
PQL<-glmmPQL(points~1,random = ~1|team,family=family,data=soccer)
Delta.start<-c(as.numeric(PQL$coef$fixed),rep(0,6),as.numeric(t(PQL$coef$random$team)))
Q.start<-as.numeric(VarCorr(PQL)[1,1])


for(j in 1:length(lambda))
{
print(paste("Iteration ", j,sep=""))
  
  for (i in 1:kk)
  {
    if (i < kk)
    {
    indi <- ind[(i-1)*nk+(1:nk)]
    }else{
    indi <- ind[((i-1)*nk+1):N]
    }
  
soccer.train<-soccer[-indi,]
soccer.test<-soccer[indi,]
  
glm2 <- try(glmmLasso(points~transfer.spendings  
        + ave.unfair.score  + ball.possession
        + tackles + ave.attend + sold.out, rnd = list(team=~1),  
        family = family, data = soccer.train, lambda=lambda[j],switch.NR=F,final.re=TRUE,
        control=list(start=Delta.start,q_start=Q.start))
        ,silent=TRUE) 
        
    if(class(glm2)!="try-error")
    {  
    y.hat<-predict(glm2,soccer.test)    

    Devianz_ma[j,i]<-sum(family$dev.resids(soccer.test$points,y.hat,wt=rep(1,length(y.hat))))
    }
}
print(sum(Devianz_ma[j,]))
}
    
Devianz_vec<-apply(Devianz_ma,1,sum)
opt2<-which.min(Devianz_vec)
       
       
glm2_final <- glmmLasso(points~transfer.spendings  
                        + ave.unfair.score + ball.possession
                        + tackles + ave.attend + sold.out, rnd = list(team=~1),  
                        family = family, data = soccer, lambda=lambda[opt2],switch.NR=F,final.re=TRUE,
                        control=list(start=Delta.start,q_start=Q.start))



summary(glm2_final)


################## More Elegant Method ############################################
## Idea: start with big lambda and use the estimates of the previous fit (BUT: before
## the final re-estimation Fisher scoring is performed!) as starting values for the next fit;
## make sure, that your lambda sequence starts at a value big enough such that all covariates are
## shrinked to zero;

## Using BIC (or AIC, respectively) to determine the optimal tuning parameter lambda
lambda <- seq(500,0,by=-5)


BIC_vec<-rep(Inf,length(lambda))
family = poisson(link = log)

# specify starting values for the very first fit; pay attention that Delta.start has suitable length! 
Delta.start<-as.matrix(t(rep(0,7+23)))
Q.start<-0.1  

for(j in 1:length(lambda))
{
  print(paste("Iteration ", j,sep=""))
  
  glm3 <- glmmLasso(points~1 +transfer.spendings
                    + ave.unfair.score 
                    + tackles  
                    + sold.out 
                    + ball.possession
                    + ave.attend
                    ,rnd = list(team=~1),  
                    family = family, data = soccer, 
                    lambda=lambda[j], switch.NR=F,final.re=TRUE,
                    control = list(start=Delta.start[j,],q_start=Q.start[j]))  
  
  print(colnames(glm3$Deltamatrix)[2:7][glm3$Deltamatrix[glm3$conv.step,2:7]!=0])
  BIC_vec[j]<-glm3$bic
  Delta.start<-rbind(Delta.start,glm3$Deltamatrix[glm3$conv.step,])
  Q.start<-c(Q.start,glm3$Q_long[[glm3$conv.step+1]])
}

opt3<-which.min(BIC_vec)

glm3_final <- glmmLasso(points~transfer.spendings + ave.unfair.score 
                        + ball.possession
                        + tackles + ave.attend + sold.out, rnd = list(team=~1),  
                        family = family, data = soccer, lambda=lambda[opt3],
                        switch.NR=F,final.re=TRUE,
                        control = list(start=Delta.start[opt3,],q_start=Q.start[opt3]))  


summary(glm3_final)

## plot coefficient paths
par(mar=c(6,6,4,4))
plot(lambda,Delta.start[2:(length(lambda)+1),2],type="l",ylim=c(-1e-1,1e-1),ylab=expression(hat(beta[j])))
lines(c(-1000,1000),c(0,0),lty=2)
for(i in 3:7){
  lines(lambda[1:length(lambda)],Delta.start[2:(length(lambda)+1),i])
}
abline(v=lambda[opt3],lty=2)


