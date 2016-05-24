library(bbmle)

set.seed(101)
z = rpois(100,lambda=5)

m1 = mle2(z~dpois(lambda=L),start=list(L=4),data=data.frame(z))

ICtab(m1,type="qAICc",dispersion=1.2,nobs=100)

m2 = glm(z~1,family=poisson)
qAICc(m2,nobs=100,dispersion=2)
