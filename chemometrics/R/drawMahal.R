"drawMahal" <-
function(x,center,covariance,quantile=c(0.975,0.75,0.5,0.25),m=1000,
   lwdcrit=1, ...){
# draw ellipses corresponding to Mahalanobis distances
# PF, 20.10.2006
#
# x ... data
# center ... center
# covariance ... covariance
# quantile: lines corresponding to quantiles of chi2
# m: number of data points on each ellipse
# lwdcrit: line width of ellipse corresponding to critical chi2 value
#

me <- center
covm <- covariance
cov.svd <- svd(covm, nv = 0)
r <- cov.svd[["u"]] %*% diag(sqrt(cov.svd[["d"]]))

alphamd <- sqrt(qchisq(quantile,2))
lalpha <- length(alphamd)

for(j in 1:lalpha) {
        e1md <- cos(c(0:m)/m * 2 * pi) * alphamd[j]
        e2md <- sin(c(0:m)/m * 2 * pi) * alphamd[j]
        emd <- cbind(e1md, e2md)
        ttmd <- t(r %*% t(emd)) + rep(1, m + 1) %o% me
        if(j == 1) {
                xmax <- max(c(x[, 1],ttmd[,1]))
                xmin <- min(c(x[, 1],ttmd[,1]))
                ymax <- max(c(x[, 2],ttmd[,2]))
                ymin <- min(c(x[, 2],ttmd[,2]))
                plot(x, xlim = c(xmin, xmax), ylim = c(ymin, ymax), ...)
        }
}
sdx <- sd(x[,1])
sdy <- sd(x[,2])
for(j in 2:lalpha) {
        e1md <- cos(c(0:m)/m * 2 * pi) * alphamd[j]
        e2md <- sin(c(0:m)/m * 2 * pi) * alphamd[j]
        emd <- cbind(e1md, e2md)
        ttmd <- t(r %*% t(emd)) + rep(1, m + 1) %o% me
        lines(ttmd[, 1], ttmd[, 2], type = "l",col=2)
}
j<-1
        e1md <- cos(c(0:m)/m * 2 * pi) * alphamd[j]
        e2md <- sin(c(0:m)/m * 2 * pi) * alphamd[j]
        emd <- cbind(e1md, e2md)
        ttmd <- t(r %*% t(emd)) + rep(1, m + 1) %o% me
        lines(ttmd[, 1], ttmd[, 2], type = "l",col=1,lwd=lwdcrit)
invisible()
}

