mnlHess =
function(beta,y,X) 
{
#   p.rossi 2004
#   changed argument order 9/05
#
# Purpose: compute mnl -Expected[Hessian]  
#
# Arguments:
#   beta is k vector of coefs
#   y is n vector with element = 1,...,j indicating which alt chosen
#   X is nj x k matrix of xvalues for each of j alt on each of n occasions
#
# Output:  -Hess evaluated at beta
#
n=length(y)
j=nrow(X)/n
k=ncol(X)
Xbeta=X%*%beta
Xbeta=matrix(Xbeta,byrow=T,ncol=j)
Xbeta=exp(Xbeta)
iota=c(rep(1,j))
denom=Xbeta%*%iota
Prob=Xbeta/as.vector(denom)

Hess=matrix(double(k*k),ncol=k)
for (i in 1:n) {
        p=as.vector(Prob[i,])
        A=diag(p)-outer(p,p)
        Xt=X[(j*(i-1)+1):(j*i),]
        Hess=Hess+crossprod(Xt,A)%*%Xt
}
return(Hess)
}
