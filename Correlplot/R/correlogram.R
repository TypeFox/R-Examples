correlogram <- function(R,labs=colnames(R),ifun="cos",cex=1,main="",ntrials=10,xlim=c(-1.2,1.2),ylim=c(-1.2,1.2),...) {
theta <- fit_angles(R,ifun=ifun,ntrials=ntrials)
X <- matrix(c(cos(theta), sin(theta)), ncol = 2)
if(is.null(labs)) labs <- 1:ncol(R)
opar <- par(pty = "s", bty="n", xaxt="n", yaxt="n")
plot(0,0, type = "n", xlab = "", ylab = "", xlim=xlim, ylim=ylim, main=main, ...)
points(X[,1],X[,2], pch = 19)
textxy(X[,1],X[,2],labs,cex=cex)
arrows(0,0,X[,1],X[,2],length=0)
circle()
par(opar)
}
