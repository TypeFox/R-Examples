library(glmm)
set.seed(1234)
data(salamander)
sal<-glmm(Mate~Cross,random=list(~0+Female,~0+Male),varcomps.names=c("F","M"), data=salamander,family.glmm=bernoulli.glmm,m=100,debug=TRUE,doPQL=FALSE)

objfun<-glmm:::objfun
beta<-rep(0,4)
nu<-rep(2,2)
par<-c(beta,nu)
nbeta<-4

debug<-sal$debug
nu.pql<-debug$nu.pql
beta.pql<-debug$beta.pql
umat<-debug$umat
m1<-debug$m1
p1<-p2<-p3<-1/3
mod.mcml<-sal$mod.mcml
u.star<-debug$u.star
family.glmm=bernoulli.glmm
zeta<-5

getEk<-glmm:::getEk
addVecs<-glmm:::addVecs

#need to get D*
eek<-getEk(mod.mcml$z)
Aks<-Map("*",eek,nu.pql)
D.star<-addVecs(Aks) 
D.star<-diag(D.star)
D.star.inv<-solve(D.star)

cache<-new.env(parent = emptyenv())

#need to also recreate the variance matrix of last imp sampling distribution
Z=do.call(cbind,mod.mcml$z)
eta.star<-as.vector(mod.mcml$x%*%beta.pql+Z%*%u.star)
cdouble<-bernoulli.glmm()$cpp(eta.star) #still a vector
cdouble<-diag(cdouble)
Sigmuh.inv<- t(Z)%*%cdouble%*%Z+D.star.inv
Sigmuh<-solve(Sigmuh.inv)

#evaluate objective function at two close places
del<-rep(10^-6,6)

ltheta<-objfun(par, nbeta, nu.pql, umat, u.star, mod.mcml, family.glmm, 
    cache, p1, p2, p3, m1, D.star, Sigmuh, Sigmuh.inv, zeta)

lthetadel<-objfun(par+del, nbeta, nu.pql, umat, u.star, mod.mcml, family.glmm, cache, p1, p2, p3, m1, D.star, Sigmuh, Sigmuh.inv, zeta) 

#do finite diffs to check value --> gradient
c(as.vector(ltheta$gradient%*%del),lthetadel$value-ltheta$value)
as.vector(ltheta$gradient%*%del)-(lthetadel$value-ltheta$value)
all.equal(as.vector(ltheta$gradient%*%del),lthetadel$value-ltheta$value,tolerance=10^-5)

#do finite diffs to check gradient --> hessian 
all.equal(lthetadel$gradient-ltheta$gradient,as.vector(ltheta$hessian%*%del),tolerance=10^-6)
cbind(lthetadel$gradient-ltheta$gradient,as.vector(ltheta$hessian%*%del))
#how big are the diffs? very small
lthetadel$gradient-ltheta$gradient-as.vector(ltheta$hessian%*%del)


