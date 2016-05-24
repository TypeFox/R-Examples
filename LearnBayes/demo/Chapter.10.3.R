#############################################################
# Section 10.3  Binary Response Regression with a Probit Link
#############################################################

#################################################
# Section 10.3.1. Missing data and Gibbs sampling
#################################################

library(LearnBayes)

 data(donner)
 attach(donner)
 X=cbind(1,age,male)

 fit=glm(survival~X-1,family=binomial(link=probit))
 summary(fit)

 m=10000
 fit=bayes.probit(survival,X,m)

 apply(fit$beta,2,mean)

 apply(fit$beta,2,sd)

 a=seq(15,65)
 X1=cbind(1,a,1)
 p.male=bprobit.probs(X1,fit$beta)

 plot(a,apply(p.male,2,quantile,.5),type="l",ylim=c(0,1),
   xlab="age",ylab="Probability of Survival")
 lines(a,apply(p.male,2,quantile,.05),lty=2)
 lines(a,apply(p.male,2,quantile,.95),lty=2)

S=readline(prompt="Type  <Return>   to continue : ")

###################################################
# Section 10.3.2  Proper priors and model selection
###################################################

library(LearnBayes)
data(donner)
y=donner$survival
X=cbind(1,donner$age,donner$male)

beta0=c(0,0,0); c0=100
P0=t(X)%*%X/c0

bayes.probit(y,X,1000,list(beta=beta0,P=P0))$log.marg

bayes.probit(y,X[,-2],1000,
   list(beta=beta0[-2],P=P0[-2,-2]))$log.marg

bayes.probit(y,X[,-3],1000,
   list(beta=beta0[-3],P=P0[-3,-3]))$log.marg

bayes.probit(y,X[,-c(2,3)],1000,
   list(beta=beta0[-c(2,3)],P=P0[-c(2,3),-c(2,3)]))$log.marg
