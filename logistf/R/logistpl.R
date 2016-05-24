logistpl <-
function(x, y, init=NULL, i, LL.0, firth, which = -1, 
offset=rep(0, length(y)), weight=rep(1,length(y)), plcontrol) {
n<-nrow(x)
k<-ncol(x)
if (is.null(init)) init<-rep(0,k)
     beta<-init
if (is.null(offset)) offset=rep(0,n)
if (is.null(weight)) weight=rep(1,n)
if (missing(plcontrol)) plcontrol<-logistpl.control()
maxit<-plcontrol$maxit
    maxstep<-plcontrol$maxstep
    maxhs<-plcontrol$maxhs
    xconv<-plcontrol$xconv
    lconv<-plcontrol$lconv
firth <- if(firth) 1 else 0
loglik <- iter <- 0
conv <- double(2)
betahist <- matrix(double(k * maxit), maxit) 
mode(x) <- mode(weight) <- mode(beta) <- mode(offset) <- mode(LL.0) <- "double"
mode(y) <- mode(firth) <- mode(n) <- mode(k) <- "integer"
mode(maxstep) <- mode(lconv) <- mode(xconv) <- mode(loglik) <- "double"
mode(maxit) <- mode(maxhs) <- mode(i) <- mode(which) <- mode(iter) <- "integer"

res <- .C("logistpl", x, y, n, k, weight, offset, beta=beta, i, which, LL.0, firth, maxit, 
maxstep, maxhs, lconv, xconv, betahist=betahist, loglik=loglik, iter=iter, conv=conv,
PACKAGE="logistf")
res <- res[c("beta", "betahist", "loglik", "iter", "conv")]
res$betahist <- head(res$betahist, res$iter)
res$beta <- res$beta[i]
res
}

