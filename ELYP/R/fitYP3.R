fitYP3 <- function(Y,d,Z,beta1=c(0,1),beta2=c(0,-1),lam=15, fun = function(t){as.numeric(t <= 0.175)})
{
#### for the given data, and given parameters
####      beta1, beta2, lam, fun compute baseline.
#### Given the baseline, compute the empirical likelihood value.
#### The number of iteration k, taken to be 16. This is the baseline iteration
#### calculation...may be a better control can be done. It seems to converge fast.

temp1 <- YP3(y=Y, d=d, Z=Z, b1=beta1, b2=beta2, k=16, lam=lam, fun=fun)
 
ELval <- ELcomp(Haz=temp1$Hazw, Sur=temp1$Survival, gam=temp1$gam)

list(LogEmpLik=ELval, MuLam=temp1$mu, BaseHazw=temp1$Hazw)
}
