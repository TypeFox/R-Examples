hiergibbs=function(data,m)
{
###############################################################
# Implements Gibbs sampling algorithm for posterior of table
# of means with hierarchical regression prior
#
# INPUT
# data:  40 by 4 matrix where the observed sample means are
#        in column 1, sample sizes are in column 2, and values of
#        two covariates in columns 3 and 4.
# m:  number of cycles of Gibbs sampling
#
# OUTPUT
# a list with 
# -- beta: matrix of simulated values of beta with each row a simulated value
# -- mu:  matrix of simulated values of cell means
# -- var:  vector of simulated values of second-stage variance sigma^2_pi
###############################################################

y=data[,1]               #
n=data[,2]               #  
x1=data[,3]              #
x2=data[,4]              # defines variables y,n,x1,x2,a 
X=cbind(1+0*x1,x1,x2)    #
s2=.65^2/n               #
p=3; N=length(y)         #

mbeta=array(0,c(m,p))        #
mmu=array(0,c(m,length(n)))  #  sets up arrays to store simulated draws
ms2pi=array(0,c(m,1))        #

########################################  defines prior parameters
b1=array(c(.55,.018,.033),c(3,1))
bvar=array(c(8.49e-03,-1.94e-05, -2.88e-04, -1.94e-05,  7.34e-07, -1.52e-06, -2.88e-04,-1.52e-06,  1.71e-05),c(3,3))
ibvar=solve(bvar)
s=.02; v=16;

mu=y; s2pi=.006  # starting values of mu and s2pi in Gibbs sampling

for (j in 1:m)
{
pvar=solve(ibvar+t(X)%*%X/s2pi)                         #
pmean=pvar%*%(ibvar%*%b1+t(X)%*%mu/s2pi)                #  simulates beta
beta=t(chol(pvar))%*%array(rnorm(p),c(p,1))+pmean       #

s2pi=(sum((mu-X%*%beta)^2)/2+s/2)/rgamma(1,shape=(N+v)/2)  #  simulates s2pi

postvar=1/(1/s2+1/s2pi)                #
postmean=(y/s2+X%*%beta/s2pi)*postvar  #  simulates mu
mu=rnorm(n,postmean,sqrt(postvar))     #

mbeta[j,]=t(beta)   #
mmu[j,]=t(mu)       #  stores simulated draws
ms2pi[j]=s2pi       #
}

return(list(beta=mbeta,mu=mmu,var=ms2pi))

}
