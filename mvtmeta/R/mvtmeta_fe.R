mvtmeta_fe <-
function(y, cov) {
if(class(y) != "matrix") 
stop("Error: input effect sizes must be put in a p by n matrix where p is the number of effects and n is the number of studies")
if(class(cov) != "array" || dim(cov)[1] != dim(cov)[2]) 
stop("Error: input covariance matrices must be put in a p by p by n array where p is the number of effects and n is the number of studies")
if(dim(y)[1] != dim(cov)[1]) 
stop("Error: the number of effects p is different in the effect size matrix and the covariance array")
if(dim(y)[2] != dim(cov)[3]) 
stop("Error: the number of studies n is different in the effect size matrix and the covariance array")
p<-dim(y)[1]
n<-dim(y)[2]
ind<-combinations(p, 2)
if(!all(sapply(1:nrow(ind), function(a) all(sapply(1:n, function(b) cov[ind[a,1],ind[a,2],b]==cov[ind[a,2],ind[a,1],b])))))
stop("Error: at least one of the covariance matrices is not symmetric")
if(!all(sapply(1:n, function(x) min(eigen(cov[,,x])$values)>0)))
stop("Error: at least one of the covariance matrices is not positive-definite")
cov_i<-array(rep(NA, p*p*n), c(p, p, n))
tmpy<-matrix(rep(NA, p*n), p, n)
for(i in 1:n) {
cov_i[,,i]<-solve(cov[,,i])
tmpy[,i]<-cov_i[,,i]%*%y[,i]
}
covbeta_i<-apply(cov_i, 1:2, sum)
covbeta<-solve(covbeta_i)
beta<-as.vector(covbeta%*%apply(tmpy, 1, sum))
return(list(beta=beta, cov=covbeta))
}

