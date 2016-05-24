"sign1" <-
function(x,makeplot=FALSE,qcrit=0.975, ...){

#################################################################
p=ncol(x)
n=nrow(x)

x.mad=apply(x,2,mad)
if (any(x.mad==0))
  stop("More than 50% equal values in one or more variables!")

# robustly sphere the data:
x.sc <- scale(x,apply(x,2,median),x.mad)
# bring to norm 1
xs <- x.sc/sqrt(apply(x.sc^2,1,sum))

# SVD; compute dimensionality:
xs.evec <- svd(xs)$v
xs.pc <- x.sc%*%xs.evec
xs.pcscal <- apply(xs.pc,2,mad)^2
xs.pcorder <- order(xs.pcscal,decreasing=TRUE)
p1=min(p-1,n-1)

# covariance matrix and distances:
covm1=xs.evec[,xs.pcorder[1:p1]]%*%diag(1/xs.pcscal[xs.pcorder[1:p1]])%*%
      t(xs.evec[,xs.pcorder[1:p1]])
x.dist=sqrt(mahalanobis(x.sc,rep(0,p),covm1,inverted=TRUE))
const <- sqrt(qchisq(qcrit,p1))
wfinal01 <- (x.dist<const)*1


#################################################################
# Generate plot:
if (makeplot){
  op <- par(mfrow=c(1,2), mar=c(4,4,2,2))
  on.exit(par(op))

  plot(x.dist,xlab="Index",ylab="Distance", ...)
  abline(h=const)
  plot(wfinal01,xlab="Index",ylab="Final 0/1 weight",ylim=c(0,1), ...)
}

list(wfinal01=wfinal01,x.dist=x.dist,const=const)
}

