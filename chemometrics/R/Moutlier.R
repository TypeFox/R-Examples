"Moutlier" <-
function(X,quantile=0.975,plot=TRUE, ...)
{
# INPUT:
# X ... data matrix
# quantile ... critical value as quantile of the chisquare distribution
# plot ... TRUE or FALSE, if TRUE plot is generated
# "..." ... additional plot parameters, see par

#require(robustbase)
X.mcd <- covMcd(X)
md=sqrt(mahalanobis(X,apply(X,2,mean),cov(X)))
rd=sqrt(mahalanobis(X,X.mcd$center,X.mcd$cov))
cutoff <- sqrt(qchisq(quantile,ncol(X)))

if (plot){
par(mfrow=c(1,2))
plot(1:length(md),md,xlab="Index of object",ylim=c(0,max(md,rd)),
  ylab="Classical Mahalanobis distance", ...)
abline(h=sqrt(qchisq(quantile,ncol(X))),lty=2)
plot(1:length(rd),rd,xlab="Index of object",ylim=c(0,max(md,rd)),
  ylab="Robust Mahalanobis distance", ...) 
abline(h=cutoff,lty=2)
}

list(md=md, rd=rd, cutoff=cutoff)
}

