fitting <- function(degree.dist,nmax){


n.regions<-length(degree.dist)
tmp<-hist(degree.dist,breaks=c(0:nmax))
cum.dist<-1-cumsum(tmp$counts)/n.regions # cumulative distribution of degree.dist

mu<-1/(sum(degree.dist)/n.regions) # parameter of the exponential law

nb<-length(degree.dist[degree.dist>0])

gamma<-1+nb/(sum(log(degree.dist[degree.dist>0]))) # parameter of the power law

x<-degree.dist
x<-x[x>0]
n<-length(x)


fn<-function(p) -(-n*p*log(sum(x)/(n*p))-n*log(gamma(p))+(p-1)*sum(log(x))-n*p)
# - log likelihood for the truncated power law with only one parameter

out <- nlm(fn, p = 1, hessian = TRUE) # Minimisation of fn

alpha<-out$estimate
beta<-sum(degree.dist)/(n.regions*alpha)
# parameters of the truncated power law




# Computation of AIC

AIC.exp<--2*(n.regions*log(mu)-mu*sum(degree.dist))+2
AIC.pow<--2*(n.regions*log(gamma-1)-gamma*sum(log(x)))+2
AIC.trunc<--2*(-out$minimum)+2

fitting<-"mu ="
fitting<-paste(fitting,mu,sep=" ")
fitting<-paste(fitting,"gamma = ",sep="\n")
fitting<-paste(fitting,gamma,sep=" ")
fitting<-paste(fitting,"alpha = ",sep="\n")
fitting<-paste(fitting,alpha,sep=" ")
fitting<-paste(fitting,"beta = ",sep="\n")
fitting<-paste(fitting,beta,sep=" ")
fitting<-paste(fitting,"AIC exp = ",sep="\n")
fitting<-paste(fitting,AIC.exp,sep=" ")
fitting<-paste(fitting,"AIC pow = ",sep="\n")
fitting<-paste(fitting,AIC.pow,sep=" ")
fitting<-paste(fitting,"AIC trunc = ",sep="\n")
fitting<-paste(fitting,AIC.trunc,sep=" ")


# file with all the value of the fittings

write.table(fitting,"fitting.txt",row.names=FALSE,col.names=FALSE,quote=FALSE)

list(mu = mu, gamma = gamma, alpha = alpha, beta = beta)

}

