##############################################
# Section 9.3 Modeling Using Zellner's g Prior
##############################################

library(LearnBayes)

# illustrating the role of the parameter c

data(puffin)
X=cbind(1, puffin$Distance - mean(puffin$Distance))
c.prior=c(0.1,0.5,5,2)
fit=vector("list",4)
for (j in 1:4)
{
  prior=list(b0=c(8,0), c0=c.prior[j])
  fit[[j]]=blinreg(puffin$Nest, X, 1000, prior)
}
BETA=NULL
for (j in 1:4)
  {
  s=data.frame(Prior=paste("c =",as.character(c.prior[j])),
         beta0=fit[[j]]$beta[,1],beta1=fit[[j]]$beta[,2])
  BETA=rbind(BETA,s)
  }
library(lattice)
with(BETA,xyplot(beta1~beta0|Prior,type=c("p","g"),col="black"))

S=readline(prompt="Type  <Return>   to continue : ")

# model selection

data=list(y=puffin$Nest, X=cbind(1,puffin$Grass,puffin$Soil))
prior=list(b0=c(0,0,0), c0=100)
beta.start=with(puffin,lm(Nest~Grass+Soil)$coef)
laplace(reg.gprior.post,c(beta.start,0),list(data=data,prior=prior))$int

X=puffin[,-1]; y=puffin$Nest; c=100
bayes.model.selection(y,X,c,constant=FALSE)
