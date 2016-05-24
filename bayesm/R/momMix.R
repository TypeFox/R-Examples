momMix=
function(probdraw,compdraw) 
{
#
# Revision History:
#   R. McCulloch 11/04
#   P. Rossi 3/05  put in backsolve fixed documentation
#   P. Rossi 9/05 fixed error in mom -- return var not sigma
#
# purpose: compute moments of normal mixture averaged over MCMC draws
#
# arguments:
#    probdraw -- ith row is ith draw of probabilities of mixture comp
#    compdraw -- list of lists of draws of mixture comp moments (each sublist is from mixgibbs)
#
# output:
#   a list with the mean vector, covar matrix, vector of std deve, and corr matrix
#
# ----------------------------------------------------------------------------------
# define function needed
mom=function(prob,comps){
# purpose: obtain mu and cov from list of normal components
#
# arguments:
#     prob: vector of mixture probs
#     comps: list, each member is a list comp with ith normal component ~N(comp[[1]],Sigma), 
#            Sigma = t(R)%*%R, R^{-1} = comp[[2]]
# returns:
#  a list with [[1]]=$mu a vector
#  [[2]]=$sigma a matrix 
#
nc = length(comps)
dim = length(comps[[1]][[1]])
mu = double(dim)
sigma = matrix(0.0,dim,dim)
for(i in 1:nc) {
   mu = mu+ prob[i]*comps[[i]][[1]]
}
var=matrix(double(dim*dim),ncol=dim)
for(i in 1:nc) {
   mui=comps[[i]][[1]]
#   root = solve(comps[[i]][[2]])
   root=backsolve(comps[[i]][[2]],diag(rep(1,dim)))
   sigma=t(root)%*%root
   var=var+prob[i]*sigma+prob[i]*(mui-mu)%o%(mui-mu)
}
list(mu=mu,sigma=var)
}
#---------------------------------------------------------------------------------------
dim=length(compdraw[[1]][[1]][[1]])
nc=length(compdraw[[1]])
dim(probdraw)=c(length(compdraw),nc)
mu=double(dim)
sigma=matrix(double(dim*dim),ncol=dim)
sd=double(dim)
corr=matrix(double(dim*dim),ncol=dim)
for(i in 1:length(compdraw)) 
{
   out=mom(probdraw[i,],compdraw[[i]])
   sd=sd+sqrt(diag(out$sigma))
   corr=corr+matrix(nmat(out$sigma),ncol=dim)
   mu=mu+out$mu
   sigma=sigma+out$sigma
}
mu=mu/length(compdraw)
sigma=sigma/length(compdraw)
sd=sd/length(compdraw)
corr=corr/length(compdraw)
return(list(mu=mu,sigma=sigma,sd=sd,corr=corr))
}
