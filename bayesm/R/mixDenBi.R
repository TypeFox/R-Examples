mixDenBi=
function(i,j,xi,xj,pvec,comps) 
{
# Revision History:
#   P. Rossi 6/05
#   vectorized evaluation of bi-variate normal density 12/06
#
# purpose: compute marg bivariate density implied by mixture of multivariate normals specified
#			by pvec,comps
#
# arguments:
#     i,j:  index of two variables
#     xi specifies a grid of points for var i
#     xj specifies a grid of points for var j
#     pvec: prior probabilities of normal components
#     comps: list, each member is a list comp with ith normal component ~ N(comp[[1]],Sigma), 
#            Sigma = t(R)%*%R, R^{-1} = comp[[2]]
# Output:
#     matrix with values of density on grid
#
# ---------------------------------------------------------------------------------------------
# define function needed
#
bivcomps=function(i,j,comps)
{
# purpose: obtain marginal means and standard deviations from list of normal components
# arguments:
#     i,j: index of elements for bivariate marginal
#     comps: list, each member is a list comp with ith normal component ~N(comp[[1]],Sigma), 
#            Sigma = t(R)%*%R, R^{-1} = comp[[2]]
# returns:
#  a list with relevant mean vectors and rooti for each compenent
#  [[2]]=$sigma a matrix whose ith row is the standard deviations for the ith component
#
result=NULL
nc = length(comps)
dim = length(comps[[1]][[1]])
ind=matrix(c(i,j,i,j,i,i,j,j),ncol=2)
for(comp in 1:nc) {
   mu = comps[[comp]][[1]][c(i,j)]
   root= backsolve(comps[[comp]][[2]],diag(dim))
   Sigma=crossprod(root)
   sigma=matrix(Sigma[ind],ncol=2)
   rooti=backsolve(chol(sigma),diag(2))
   result[[comp]]=list(mu=mu,rooti=rooti)
}
return(result)
}
# ----------------------------------------------------------------------------------------------
nc = length(comps)
marmoms=bivcomps(i,j,comps)
ngridxi=length(xi); ngridxj=length(xj)
z=cbind(rep(xi,ngridxj),rep(xj,each=ngridxi))
den = matrix(0.0,nrow=ngridxi,ncol=ngridxj)
for(comp in 1:nc) {
  quads=colSums((crossprod(marmoms[[comp]]$rooti,(t(z)-marmoms[[comp]]$mu)))^2)
  dencomp=exp(-(2/2)*log(2*pi)+sum(log(diag(marmoms[[comp]]$rooti)))-.5*quads) 
  dim(dencomp)=c(ngridxi,ngridxj)
  den=den+dencomp*pvec[comp]
}
return(den)
}
