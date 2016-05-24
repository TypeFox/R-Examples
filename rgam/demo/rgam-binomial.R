
library(rgam)
library(gam)
library(akima)

# Generate a data set
n <- 2000
set.seed(123456)
x1 <- rnorm(n, mean=1/2, sd=.2)
x2 <- rexp(n, rate=.6)
nu <- sin(pi*x1)*5 - log(x2+2)*4
pr <- 1 / (1+exp(-nu))
y <- rbinom(n, size=10, prob=pr)
ind <- sample(1:n, 300)
x1 <- x1[ind]
x2 <- x2[ind]
y <- y[ind]
pr <- pr[ind]
n <- length(x1)
ni <- rep(10, n)

### add outliers
i <- c(14,  35,  89, 98, 106, 115, 158, 181, 205, 290, 10, 102, 267, 283)
y[i] <- floor(ni[i]*.9)


### Robust and standard fits
a <- rgam(x=cbind(x1, x2), y=y, ni=ni, family='binomial', cv.method='rcv', alpha=.3, trace=TRUE)
b <- gam(cbind(y, ni-y)~lo(x1, span=.3) + lo(x2, span=.3), family=binomial)

### Predicted probabilities of success
pr.gam <- as.vector(predict(b, type='response'))
pr.rgam <- predict(a, type='response')

### interpolate the predicted points (for display purposes only)
no <- 100
dd <- interp(x=x1, y=x2, z=pr.gam, xo=seq(min(x1), max(x1), length = no),
yo=seq(min(x2), max(x2), length = no),
linear=TRUE, duplicate='mean')

### display the predicted surfaces and the data
par(mfrow=c(1,2))
pp <- persp(dd$x, dd$y, matrix(dd$z, no, no), theta=120, phi=10, 
zlim=c(0, 1), col='lightblue', shade=.75, main='GAM fit', sub='gray points are below the surface')
vv <- ((pr.gam - y/ni)>0)
points(trans3d(x1, x2, (y/ni), pmat=pp), pch=19, cex=.75, col=c('red', 'gray90')[vv+1])

no <- 100
dd.r <- interp(x=x1, y=x2, z=pr.rgam, xo=seq(min(x1), max(x1), length = no),
yo=seq(min(x2), max(x2), length = no),
linear=TRUE, duplicate='mean')

pp.r <- persp(dd.r$x, dd.r$y, matrix(dd.r$z, no, no), theta=120, phi=10, 
zlim=c(0, 1), col='lightblue', shade=.75, main='Robust GAM fit', sub='gray points are below the surface')

vv <- ((pr.rgam - y/ni)>0)
points(trans3d(x1, x2, (y/ni), pmat=pp.r), pch=19, cex=.75, col=c('red', 'lightgray')[vv+1])

par(mfrow=c(1,1))


### residual plots
vv <- (pr.gam - y/ni) / sqrt(pr.gam*(1-pr.gam))
vv.r <- (pr.rgam - y/ni) / sqrt(pr.rgam*(1-pr.rgam))
par(mfrow=c(1,2))
plot(vv, ylim=c(-7, 1.5), pch=19, col='gray', ylab='Standardized residual', main='GAM fit')
abline(h=0, lwd=4, col='lightblue')
plot(vv.r, ylim=c(-7, 1.5), pch=19, col='gray', ylab='Standardized residual', main='Robust GAM fit')
abline(h=0, lwd=4, col='lightblue')
par(mfrow=c(1,1))

### true versus predicted plots
par(mfrow=c(1,2))
plot(pr, pr.gam, pch=19, col='gray', xlab='True probabilities', ylab='Predicted values', main='GAM fit', xlim=c(0,1), ylim=c(0,1))
abline(0,1, lwd=4, col='lightblue')
plot(pr, pr.rgam, pch=19, col='gray', xlab='True probabilities', ylab='Predicted values', main='Robust GAM fit', xlim=c(0,1), ylim=c(0,1))
abline(0,1, lwd=4, col='lightblue')
par(mfrow=c(1,1))


