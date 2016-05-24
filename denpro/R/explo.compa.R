explo.compa<-function(dendat,seed=1)
{
n<-dim(dendat)[1]
d<-dim(dendat)[2]

cova<-cov(dendat)
mu<-mean(data.frame(dendat))

eig<-eigen(cov(dendat),symmetric=TRUE)
sigsqm<-eig$vectors%*%diag(eig$values^{1/2})  #%*%t(eig$vectors)

set.seed(seed)
symmedata<-matrix(rnorm(d*n),n,d)
dendat.simu<-t(sigsqm%*%t(symmedata))

return(dendat.simu)
}

