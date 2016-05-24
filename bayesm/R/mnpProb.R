mnpProb=
function(beta,Sigma,X,r=100) 
{
#
# revision history:
#  written by Rossi 9/05
#  W. Taylor 4/15 - replaced ghkvec call with rcpp version
#
# purpose:
#   function to MNP probabilities for a given X matrix (corresponding
#   to "one" observation
#
# arguments: 
#   X is p-1 x k array of covariates (including intercepts)
#      note: X is from the "differenced" system
#   beta is k x 1  with k = ncol(X)
#   Sigma is p-1 x p-1 
#   r is the number of random draws to use in GHK
#
# output -- probabilities
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
pm1=ncol(Sigma)
k=length(beta)
mu=matrix(X%*%beta,nrow=pm1)
above=rep(0,pm1)
prob=double(pm1+1)
for (j in 1:pm1) {
   Aj=-diag(pm1)
   Aj[,j]=rep(1,pm1)
   trunpt=as.vector(-Aj%*%mu)
   Lj=t(chol(Aj%*%Sigma%*%t(Aj)))
#     note: rob's routine expects lower triangular root
   prob[j]=ghkvec(Lj,trunpt,above,r)
#     note:  ghkvec does an entire vector of n probs each with different truncation points but the
#            same cov matrix.  
}
#
# now do pth alternative
#
prob[pm1+1]=1-sum(prob[1:pm1])
return(prob)

}
