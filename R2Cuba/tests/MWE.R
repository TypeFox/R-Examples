# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Example 2015/09/03
# Pierre de Villemereuil : pierre.de.villemereuil@mailoo.org
# Call of "cuhre".
# The integrand calls itself "cuhre"
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


library(R2Cuba)
library(mvtnorm)

#Minimal working example
mu=c(0,10)			#Mean
G=matrix(c(0.5,0,0,1),nrow=2)	#Some variance-covariance matrix
P=matrix(c(1,0,0,2),nrow=2) 	#Some other VCV matrix

#Arbitrary function yielding a scalar
arb.func<-function(x){x[1]+0.5*x[2]}
rel.tol<-10; abs.tol<-10# Not so good, for timing purpose


#We want to compute the covariance between a vector v and an arbitrary function of another vector which depends on v
#A way to do that is first to compute the expectency of the function given a value of the vector v
exp_func_v<-function(v){cuhre(ndim=2,ncomp=1,integrand=function(x){arb.func(x)*dmvnorm(x,mu+v,P)},lower=-c(100,1000),upper=c(100,1000),
rel.tol=rel.tol, abs.tol=abs.tol,
flags=list(verbose=0))$value}
#exp_func_v works and indeed yields a scalar
exp_func_v(c(10,0))

#Then we average the expectancy above over all values of v
a<-cuhre(ndim=2,ncomp=2,integrand=function(v){v*exp_func_v(v)*dmvnorm(v,c(0,0),G)}, lower=-c(100,1000),upper=c(100,1000),
rel.tol=rel.tol, abs.tol=abs.tol,
flags=list(verbose=3))
print(a$value)

