myLLfun <- function(mle, dataMat, fun=function(t){as.numeric(t <= 0.25)}) {
#### mle is a vector of length 3. the parameters beta1 beta2 and lam.
####   Technically, the mle of lam is (always) 0. But in this function we allow 
####    other parameter values as input; i.e. the input beta1, beta2 lam do not have to be mle.
####   This function also returns Mulan, which is a 1-1 function of lam. The definition of Mulam
####    is controlled by the fun input in the fitYP3( ). 
#### dataMat  should be a matrix of 3 by n, which is Y,d,Z vectors. ???? or 4 by n but fourth row do not matter.
#### It returns the MuLam = fun(lam), and the log EmpLik value.
b1 <- mle[1]
b2 <- mle[2]
lam <- mle[3]
ftemp <- fitYP3(Y=dataMat[1,], d=dataMat[2,], Z=t(dataMat[3:4,]), beta1=c(0,b1), beta2=c(0,b2), lam=lam, fun=fun)
list(Mulam=ftemp$MuLam, Loglik=ftemp$LogEmpLik)
}
