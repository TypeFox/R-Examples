myLLfun2 <- function(mle, dataMat, fun=function(t){as.numeric(t <= 0.25)}) {
#### mle is a vector of length 4. the parameters beta1 beta2 a and lam.
####   Technically, the mle of lam is (always) 0. But in this function we allow 
####    other parameter values as input; i.e. the input beta1, beta2, a, lam do not have to be mle.
####   This function also returns Mulan, which is a 1-1 function of lam. The definition of Mulam
####    is controlled by the fun input. 
#### dataMat  should be a matrix of 4 by n, which is Y,d,Z vectors.
#### It returns the MuLam = fun(lam), and the log EmpLik value.
#### The difference to myLLfun( ) is that we have one more input: a, and the dataMat also one more row.
####  This function is used by findU4( ), findL4
####
b1 <- mle[1]
b2 <- mle[2]
a <- mle[3]
lam <- mle[4]
ftemp <- fitYP3(Y=dataMat[1,], d=dataMat[2,], Z=t(dataMat[3:4,]), beta1=c(a,b1), beta2=c(a,b2), lam=lam, fun=fun)
list(Mulam=ftemp$MuLam, Loglik=ftemp$LogEmpLik)
}
