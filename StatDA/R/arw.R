"arw" <-
function(x,m0,c0,alpha,pcrit){
# Adaptive reweighted estimator for multivariate location and scatter
# with hard-rejection weights and delta = chi2inv(1-d,p)
#
# Input arguments
#   x:  Dataset (n x p)
#   m0: Initial location estimator (1 x p)
#   c0: Initial scatter estimator (p x p)
#   alpha:  Maximum thresholding proportion
#       (optional scalar, default: alpha = 0.025)
#   pcrit: critical value for outlier probability
#       (optional scalar, default values from simulations)
#
# Output arguments:
#   m:  Adaptive location estimator (p x 1)
#   c:  Adaptive scatter estimator (p x p)
#   cn: Adaptive threshold (scalar)
#   w:  Weight vector (n x 1)
#
n <- nrow(x)
p <- ncol(x)
# Critical value for outlier probability based on simulations for alpha=0.025
if (missing(pcrit)){
  if (p<=10) pcrit <- (0.24-0.003*p)/sqrt(n)
  if (p>10) pcrit <- (0.252-0.0018*p)/sqrt(n)
}
if (missing(alpha)) delta<-qchisq(0.975,p) else delta<-qchisq(1-alpha,p)
d2<-mahalanobis(x,m0,c0)
d2ord <- sort(d2)
dif <- pchisq(d2ord,p) - (0.5:n)/n
i <- (d2ord>=delta) & (dif>0)
if (sum(i)==0) alfan<-0 else alfan<-max(dif[i])
if (alfan<pcrit) alfan<-0
#if (alfan>0) cn<-max(d2ord[n-floor(n*alfan)],delta) else cn<-Inf
if (alfan>0) cn<-max(d2ord[n-ceiling(n*alfan)],delta) else cn<-Inf
w <- d2<cn
if(sum(w)==0) {
  m <- m0
  c <- c0
} 
else {
  m <- apply(x[w,],2,mean)
  c1 <- as.matrix(x-rep(1,n)%*%t(m))
  c <- (t(c1)%*%diag(w)%*%c1)/sum(w)
}
list(m=m,c=c,cn=cn,w=w)
}

