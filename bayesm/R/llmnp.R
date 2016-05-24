llmnp=
function(beta,Sigma,X,y,r) 
{
#
# revision history:
#   edited by rossi 2/8/05
#   adde 1.0e-50 before taking log to avoid -Inf 6/05
#   changed order of arguments to put beta first 9/05
#
# purpose:
#   function to evaluate MNP likelihood using GHK
#
# arguments: 
#   X is n*(p-1) x k array of covariates (including intercepts)
#      note: X is from the "differenced" system
#   y is vector of n indicators of multinomial response
#   beta is k x 1  with k = ncol(X)
#   Sigma is p-1 x p-1 
#   r is the number of random draws to use in GHK
#
# output -- value of log-likelihood
# for each observation w = Xbeta + e   e ~N(0,Sigma)
#   if y=j (j<p)
#      w_j > max(w_-j) and w_j >0
#   if y=p,  w < 0
#
# to use GHK we must transform so that these are rectangular regions
# e.g.  if y=1,  w_1 > 0 and w_1 - w_-1 > 0 
#
# define Aj such that if j=1,..,p-1, Ajw = Ajmu + Aje > 0 is equivalent to y=j
#   implies Aje > -Ajmu
#   lower truncation is -Ajmu and cov = AjSigma t(Aj)
#
# for p, e < - mu
#
#
# define functions needed
#
ghkvec = function(L,trunpt,above,r){
   dim=length(above)
   n=length(trunpt)/dim
   .C('ghk_vec',as.integer(n),as.double(L),as.double(trunpt),as.integer(above),as.integer(dim),
   as.integer(r),res=double(n))$res}
#   
# compute means for each observation
#
pm1=ncol(Sigma)
k=length(beta)
mu=matrix(X%*%beta,nrow=pm1)
logl=0.0
above=rep(0,pm1)
for (j in 1:pm1) {
   muj=mu[,y==j]
   Aj=-diag(pm1)
   Aj[,j]=rep(1,pm1)
   trunpt=as.vector(-Aj%*%muj)
   Lj=t(chol(Aj%*%Sigma%*%t(Aj)))
#     note: rob's routine expects lower triangular root
   logl=logl + sum(log(ghkvec(Lj,trunpt,above,r)+1.0e-50))
#     note:  ghkvec does an entire vector of n probs each with different truncation points but the
#            same cov matrix.  
}
#
# now do obs for y=p
#
trunpt=as.vector(-mu[,y==(pm1+1)])
Lj=t(chol(Sigma))
above=rep(1,pm1)
logl=logl+sum(log(ghkvec(Lj,trunpt,above,r)+1.0e-50))
return(logl)

}
