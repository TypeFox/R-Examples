mixDen=
function(x,pvec,comps) 
{
# Revision History:
#   R. McCulloch 11/04
#   P. Rossi 3/05 -- put in backsolve
#   P. Rossi 1/06 -- put in crossprod
#
# purpose: compute marginal densities for multivariate mixture of normals (given by p and comps) at x
#
# arguments:
#     x: ith columns gives evaluations for density of ith variable
#     pvec: prior probabilities of normal components
#     comps: list, each member is a list comp with ith normal component ~ N(comp[[1]],Sigma), 
#            Sigma = t(R)%*%R, R^{-1} = comp[[2]]
# Output:
#     matrix with same shape as input, x, ith column gives margial density of ith variable
#
# ---------------------------------------------------------------------------------------------
# define function needed
#
ums=function(comps)
{
# purpose: obtain marginal means and standard deviations from list of normal components
# arguments:
#     comps: list, each member is a list comp with ith normal component ~N(comp[[1]],Sigma), 
#            Sigma = t(R)%*%R, R^{-1} = comp[[2]]
# returns:
#  a list with [[1]]=$mu a matrix whose ith row is the means for ith component
#  [[2]]=$sigma a matrix whose ith row is the standard deviations for the ith component
#
nc = length(comps)
dim = length(comps[[1]][[1]])
mu = matrix(0.0,nc,dim)
sigma = matrix(0.0,nc,dim)
for(i in 1:nc) {
   mu[i,] = comps[[i]][[1]]
#   root = solve(comps[[i]][[2]])
   root= backsolve(comps[[i]][[2]],diag(rep(1,dim)))
   sigma[i,] = sqrt(diag(crossprod(root)))
}
return(list(mu=mu,sigma=sigma))
}
# ----------------------------------------------------------------------------------------------
nc = length(comps)
mars = ums(comps)
den = matrix(0.0,nrow(x),ncol(x))
for(i in 1:ncol(x)) {
   for(j in 1:nc) den[,i] = den[,i] + dnorm(x[,i],mean = mars$mu[j,i],sd=mars$sigma[j,i])*pvec[j]
}
return(den)
}
