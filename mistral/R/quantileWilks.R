quantileWilks = function(alpha=0.95, beta=0.95, data=NULL, bilateral=FALSE){
  
# Input:     2 required arguments (by default: unilateral quantile)
# alpha:     unilateral or bilateral quantile level
# beta:      confidence level on quantile value(s)
# data:      the data sample (vector format) to compute the quantile(s)
#            if data=NULL, the function returns the minimal sample size to compute the required quantile 
# bilateral: TRUE  -> bilateral quantiles (prediction/tolerance interval) of order alpha
#            FALSE -> unilateral quantiles (right size) of order alpha
#
# Output: 4 values if 'data' is specified
#         1 value (nmin) if 'data' is not specified
# lower:  lower bound of the bilateral tolerance interval
#         if bilateral=FALSE, no value 
# upper:  upper bound of the tolerance interval (bilateral case) or quantile value  (unilateral case)
# nmin:   minimal size of the required i.i.d. sample for given alpha and beta 
#           - bilateral case: tolerance interval will be composed with the min and max of the sample
#           - unilateral case: the quantile will correspond to max of the sample
# ind:    the index (unilateral case) or indices (bilateral case) of the quantiles 
#         in the ordered sample (increasing order) 

if (!bilateral){
f = function(n){1-alpha^n-beta}
nmin = ceiling(uniroot(f,c(0,1e7))$root)
}

if (bilateral){
f = function(n){1-alpha^n-n*(1-alpha)*alpha^(n-1)-beta}
nmin = ceiling(uniroot(f,c(0,1e7))$root)
}

n = length(data)

if(n!=0){

if (n<nmin){
cat('n < nmin, change of alpha and/or beta values \n')
return(list(nmin=nmin))
}

if (n>=nmin){

if (bilateral){
up = floor(n/2)
rseq = 1:up
out = cbind(rseq,pbeta(alpha,n-2*rseq+1,2*rseq,lower.tail=FALSE) - beta)
r = which.min(out[which(out[,2]>0),2])
upper = sort(data,decreasing=TRUE)[r]
lower = sort(data,decreasing=FALSE)[r]
output = list(lower=lower,upper=upper,nmin=nmin,ind=c(r,n-r+1))
return(output)
}

if (!bilateral){
r = qbinom(1-beta,n,1-alpha)
upper = sort(data,decreasing=FALSE)[n-r+1]
output = list(upper=upper,nmin=nmin,ind=n-r+1)
return(output)
}
}
}
else{return(list(nmin=nmin))}
}
