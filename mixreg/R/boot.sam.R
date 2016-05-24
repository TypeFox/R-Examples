boot.sam <- function(mu,resid,prob) {
#
# Function boot.sam to draw a bootstrap sample from the (multiple)
# residuals of a mixture of regressions model.
#
# Think of the model as follows:  Marbles are labelled 1, 2, ..., n;
# each marble contains a set of pairs (j, eps_j) for j = 1, 2, ..., k
# where k is the number of components.  Each pair is labelled with a
# probability p_j, with p_1 + ... + p_k = 1 of course.  To form y_i^*
# choose a marble with probability 1/n; then choose a pair with
# probability p_j.  Then set y_i^* = mu_ij + eps_j.
#

nc  <- ncol(mu)
n   <- nrow(mu)
if(nc==1) return(drop(mu+sample(resid,n,TRUE)))

ii  <- sample(1:n,n,TRUE)
jj  <- apply(prob[ii,],1,function(x,k){sample(1:k,1,prob=x)},nc)
drop(mu[cbind(1:n,jj)] + resid[cbind(ii,jj)])

}
