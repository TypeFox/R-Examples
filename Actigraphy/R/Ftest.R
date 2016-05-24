Ftest <-
function(yfdPar, xfdlist, betalist, wt=NULL, nperm=1000, argvals=NULL, q=0.05, plotres=TRUE, mul=1, ...){
if(missing(yfdPar) || missing(xfdlist) || missing(betalist)) 
stop("Missing Arguments")

set.seed <- 123456789
Fnull <- rep(0, nperm)
Fnullvals <- c()
q <- 1 - q
begin <- proc.time()
fRegressList <- fRegress(yfdPar, xfdlist, betalist)
elapsed.time <- max(proc.time() - begin, na.rm=TRUE)     
print(paste("Permutation F test running (", nperm, " permutations)", sep=""))
print(paste("Estimated Computing time = ", round(nperm * elapsed.time), " seconds", sep=""))

yhat <- fRegressList$yhatfdobj
if(is.list(yhat) && ("fd" %in% names(yhat))) 
yhat <- yhat$fd
tFstat <- Fstat.fd(yfdPar, yhat, argvals)
Fvals <- tFstat$F
Fobs <- max(Fvals)
argvals <- tFstat$argvals

if(is.vector(yfdPar)){
n <- length(yfdPar)
}else{
n <- ncol(yfdPar$coefs)
}

for(i in 1:nperm){
tyfdPar <- yfdPar[sample(n)]
fRegressList <- fRegress(tyfdPar, xfdlist, betalist)
yhat <- fRegressList$yhatfdobj
if(is.list(yhat) && ("fd" %in% names(yhat))) 
yhat <- yhat$fd
tFstat <- Fstat.fd(yfdPar, yhat, argvals)
Fnullvals <- cbind(Fnullvals, tFstat$F)
Fnull[i] <- max(Fnullvals[, i])
}

pval <- mean(Fobs < Fnull)
qval <- quantile(Fnull, q)
pvals.pts <- apply(Fvals < Fnullvals, 1, mean)
qvals.pts <- apply(Fnullvals, 1, quantile, q)
if(plotres){
if(is.fd(yfdPar)){
ylims <- c(min(c(Fvals, qval, qvals.pts)), max(c(Fobs, qval))*mul)

if(is.null(names(yhat$fdnames))){
xlab <- "argvals"
}else{
xlab <- names(yhat$fdnames)[1]
}

plot(argvals, Fvals, type="l", ylim=ylims, col=2, lwd=2, xlab=" ", ylab="F-statistic", main="Permutation F-Test", ...)
lines(argvals, qvals.pts, lty=3, col=4, lwd=2)
abline(h=qval, lty=2, col=4, lwd=2)
legendstr <- c("Observed Statistic", paste("pointwise", 1 - q, "critical value"), paste("maximum", 1 - q, "critical value"))
legend("topleft", legend=legendstr, col=c(2,4,4), lty=c(1, 3, 2), lwd=c(2,2,2), cex=0.8)
}else{
xlims <- c(min(c(Fnull, Fobs)), max(c(Fnull, Fobs)))
hstat <- hist(Fnull, xlim=xlims, lwd=2, xlab="F-value", main="Permutation F-Test")
abline(v=Fobs, col=1, lwd=2)
abline(v=qval, col=1, lty=2, lwd=2)
legendstr <- c("Observed Statistic", paste("Permutation", 1 - q, "critical value"))
legend("topleft", legend=legendstr, lty=c(1, 2), lwd=c(2,2))
}
}

return(list(pval=pval, qval=qval, Fobs=Fobs, Fnull=Fnull, Fvals=Fvals, 
Fnullvals=Fnullvals, pvals.pts=pvals.pts, qvals.pts=qvals.pts, 
fRegressList=fRegressList, argvals=argvals))
}
