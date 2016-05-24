### R code from vignette source 'profdpm.Rnw'

###################################################
### code chunk number 1: profdpm.Rnw:194-195
###################################################
options(prompt = "R> ", continue = "+  ", width = 70, useFancyQuotes = FALSE)


###################################################
### code chunk number 2: profdpm.Rnw:197-212
###################################################
set.seed(42)

sim <- function(multiplier = 1) {
    x <- as.matrix(runif(99))
    a <- multiplier * c(5,0,-5)
    s <- multiplier * c(-10,0,10)
    y <- c(a[1]+s[1]*x[1:33],
           a[2]+s[2]*x[34:66],
           a[3]+s[3]*x[67:99]) + rnorm(99)
    group <- rep(1:33, rep(3,33))
    return(data.frame(x=x,y=y,gr=group))
}
dat <- sim()
library("profdpm")
fitL <- profLinear(y ~ x, group=gr, data=dat)


###################################################
### code chunk number 3: figSim1 (eval = FALSE)
###################################################
## sfitL <- summary(fitL)
## plot(fitL$x[,2], fitL$y, col=grey(0.9), xlab='x', ylab='y')
## for(grp in unique(fitL$group)) {
##     ind <- which(fitL$group==grp)
##     ord <- order(fitL$x[ind,2])
##     lines(fitL$x[ind,2][ord],
##           fitL$y[ind][ord],
##           col=grey(0.9))
## }
## for(cls in 1:length(sfitL)) {
##     # The following implements the (3rd) method of
##     # Hanson & McMillan (2012) for simultaneous credible bands
##     
##     # Generate coefficients from profile posterior
##     n   <- 1e4
##     tau <- rgamma(n, shape=fitL$a[[cls]]/2, scale=2/fitL$b[[cls]])
##     muz <- matrix(rnorm(n*2, 0, 1),n,2)
##     mus <- (muz / sqrt(tau)) %*% chol(solve(fitL$s[[cls]]))
##     mu  <- outer(rep(1,n), fitL$m[[cls]]) + mus  
##     
##     # Compute Mahalanobis distances
##     mhd <- rowSums(muz^2)
##     
##     # Find the smallest 95% in terms of Mahalanobis distance
##     # I.e., a 95% credible region for mu
##     ord <- order(mhd, decreasing=TRUE)[-(1:floor(n*0.05))]
##     mu  <- mu[ord,]
##     
##     #Compute the 95% credible band
##     plotx <- seq(min(dat$x), max(dat$x), length.out=200)
##     ral  <- apply(mu, 1, function(m) m[1] + m[2] * plotx)
##     rlo  <- apply(ral, 1, min)
##     rhi  <- apply(ral, 1, max)
##     rmd  <- fitL$m[[cls]][1] + fitL$m[[cls]][2] * plotx
##     
##     lines(plotx, rmd, col=cls, lty=2)
##     lines(plotx, rhi, col=cls)
##     lines(plotx, rlo, col=cls)
## }


###################################################
### code chunk number 4: profdpm.Rnw:262-263
###################################################
sfitL <- summary(fitL)
plot(fitL$x[,2], fitL$y, col=grey(0.9), xlab='x', ylab='y')
for(grp in unique(fitL$group)) {
    ind <- which(fitL$group==grp)
    ord <- order(fitL$x[ind,2])
    lines(fitL$x[ind,2][ord],
          fitL$y[ind][ord],
          col=grey(0.9))
}
for(cls in 1:length(sfitL)) {
    # The following implements the (3rd) method of
    # Hanson & McMillan (2012) for simultaneous credible bands
    
    # Generate coefficients from profile posterior
    n   <- 1e4
    tau <- rgamma(n, shape=fitL$a[[cls]]/2, scale=2/fitL$b[[cls]])
    muz <- matrix(rnorm(n*2, 0, 1),n,2)
    mus <- (muz / sqrt(tau)) %*% chol(solve(fitL$s[[cls]]))
    mu  <- outer(rep(1,n), fitL$m[[cls]]) + mus  
    
    # Compute Mahalanobis distances
    mhd <- rowSums(muz^2)
    
    # Find the smallest 95% in terms of Mahalanobis distance
    # I.e., a 95% credible region for mu
    ord <- order(mhd, decreasing=TRUE)[-(1:floor(n*0.05))]
    mu  <- mu[ord,]
    
    #Compute the 95% credible band
    plotx <- seq(min(dat$x), max(dat$x), length.out=200)
    ral  <- apply(mu, 1, function(m) m[1] + m[2] * plotx)
    rlo  <- apply(ral, 1, min)
    rhi  <- apply(ral, 1, max)
    rmd  <- fitL$m[[cls]][1] + fitL$m[[cls]][2] * plotx
    
    lines(plotx, rmd, col=cls, lty=2)
    lines(plotx, rhi, col=cls)
    lines(plotx, rlo, col=cls)
}


###################################################
### code chunk number 5: profdpm.Rnw:278-281
###################################################
simulatedPartL <- rep(1:3, rep(33,3))
estimatedPartL <- fitL$clust
pci(simulatedPartL, estimatedPartL)


###################################################
### code chunk number 6: figSim2 (eval = FALSE)
###################################################
## mult <- rep(seq(1, 0.01, length.out=33), 10)
## Rand <- sapply(mult, function(mult) {
##     dat <- sim(mult)
##     fit <- profLinear(y ~ x, group=gr, data=dat) 
##     pci(rep(1:3, rep(33,3)), fit$clust)['R']
## })
## lws <- lowess(x=mult,y=Rand)
## plot(lws$x, lws$y, type='n',
##      xlab='coefficient multiplier', ylab='Rand index')
## points(mult + rnorm(length(mult),0,1/200),
##      Rand + rnorm(length(Rand),0,1/200), col=grey(0.9))
## lines(lws$x, lws$y)


###################################################
### code chunk number 7: profdpm.Rnw:309-310
###################################################
mult <- rep(seq(1, 0.01, length.out=33), 10)
Rand <- sapply(mult, function(mult) {
    dat <- sim(mult)
    fit <- profLinear(y ~ x, group=gr, data=dat) 
    pci(rep(1:3, rep(33,3)), fit$clust)['R']
})
lws <- lowess(x=mult,y=Rand)
plot(lws$x, lws$y, type='n',
     xlab='coefficient multiplier', ylab='Rand index')
points(mult + rnorm(length(mult),0,1/200),
     Rand + rnorm(length(Rand),0,1/200), col=grey(0.9))
lines(lws$x, lws$y)


###################################################
### code chunk number 8: profdpm.Rnw:329-334
###################################################
p <- seq(0.9,0.1,length.out=9)
y1 <- matrix(rbinom(999, 1, p), 111, 9, TRUE)
y2 <- matrix(rbinom(999, 1, rev(p)), 111, 9, TRUE)
dat <- as.data.frame(rbind(y1, y2))
fitb <- profBinary(~0+., data=dat)


###################################################
### code chunk number 9: figSim3
###################################################
alpha <- 10^(-seq(0, 40, length.out=33))
Rand <- sapply(alpha, function(alpha) {
    fit <- profBinary(~0+., data=dat, param=list(alpha=alpha))
    pci(rep(1:2, rep(111,2)), fit$clust)['R']
})
lws <- lowess(x=log(alpha),y=Rand)
plot(lws$x, lws$y, type='l',
     ylim=c(0.45,1),
     xlab=expression(paste('log', alpha)),
     ylab='Rand index')

points(log(alpha), Rand + rnorm(length(Rand),0,1/200),
     col=grey(0.9))


###################################################
### code chunk number 10: profdpm.Rnw:362-363
###################################################
alpha <- 10^(-seq(0, 40, length.out=33))
Rand <- sapply(alpha, function(alpha) {
    fit <- profBinary(~0+., data=dat, param=list(alpha=alpha))
    pci(rep(1:2, rep(111,2)), fit$clust)['R']
})
lws <- lowess(x=log(alpha),y=Rand)
plot(lws$x, lws$y, type='l',
     ylim=c(0.45,1),
     xlab=expression(paste('log', alpha)),
     ylab='Rand index')

points(log(alpha), Rand + rnorm(length(Rand),0,1/200),
     col=grey(0.9))


