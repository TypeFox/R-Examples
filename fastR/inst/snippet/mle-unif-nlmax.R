x <- c(1.6,2.8,6.2,8.2,8.7)
loglik <- function(theta,x) {
     sum ( dunif(x,0,theta,log=TRUE) )
}
lik <- function(theta,x) {
     prod ( dunif(x,0,theta,log=FALSE) )
}

# works just fine if we select a good starting point.
summary(nlmax(loglik, p=10, x=x))
# but some starting points don't work well...
p <- seq(7,12,by=0.2)
est <- sapply(p, function(p) { nlmax(loglik,p=p, x=x)$estimate })
rbind(p,est)
# here's another try without the logarithmic transformation.
est <- sapply(p, function(p) { nlmax(lik,p=p, x=x)$estimate })
rbind(p,est)
# a graph of the likelihood function shows why
theta <- seq(6,12,by=0.002)
y1 <- sapply(theta, function(theta) { lik(theta,x)} )
y2 <- sapply(theta, function(theta) { loglik(theta,x)} )
plot1 <- xyplot(y1 ~ theta,
    xlab=expression(theta),
    ylab='likelihood',cex=0.5)
plot2 <- xyplot(y2 ~ theta,
    xlab=expression(theta),
    ylab='log-likelihood',cex=0.5)
