"sign2" <-
function(x,makeplot=FALSE,explvar=0.99,qcrit=0.975, ...){

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

# PC decomposition; compute p*, robustly sphere:
xs.evec <- svd(xs)$v
xs.pc <- x.sc%*%xs.evec
xs.pcscal <- apply(xs.pc,2,mad)^2
xs.pcorder <- order(xs.pcscal,decreasing=TRUE)
p1 <- (1:p)[(cumsum(xs.pcscal[xs.pcorder])/sum(xs.pcscal)>explvar)][1]
x.pc <- x.sc%*%xs.evec[,xs.pcorder[1:p1]]
xpc.sc <- scale(x.pc,apply(x.pc,2,median),apply(x.pc,2,mad))

# compute distances and weights
xpc.norm <- sqrt(apply(xpc.sc^2,1,sum))
xpc.out <- xpc.norm/median(xpc.norm)
x.dist <- xpc.out*sqrt(qchisq(0.5,p1))
const <- sqrt(qchisq(qcrit,p1))
wfinal01 <- rep(0,n)
wfinal01[x.dist<const] <- 1


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

