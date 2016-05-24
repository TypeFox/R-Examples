### R code from vignette source 'article.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: initial_settings
###################################################
options(stringsAsFactors=FALSE)
options(width=65)
options(prompt="> ")


###################################################
### code chunk number 2: eq1a
###################################################
library("genlasso")
set.seed(1)
n = 100
p = 10
X = matrix(rnorm(n*p), ncol=p)
y = X[,1] + rnorm(n)


###################################################
### code chunk number 3: eq1b
###################################################
D = diag(1,p)


###################################################
### code chunk number 4: eq1c
###################################################
out = genlasso(y, X=X, D=D)


###################################################
### code chunk number 5: eq1d
###################################################
out


###################################################
### code chunk number 6: eq1e
###################################################
summary(out)


###################################################
### code chunk number 7: eq1f
###################################################
plot(out)


###################################################
### code chunk number 8: figeq1f
###################################################
plot(out)


###################################################
### code chunk number 9: article.Rnw:251-252
###################################################
coef(out, lambda=sqrt(n*log(p)))


###################################################
### code chunk number 10: eq2a
###################################################
set.seed(1)
n = 100
i = 1:n
y = (i > 20 & i < 30) + 5*(i > 50 & i < 70) +
  rnorm(n, sd=0.1)


###################################################
### code chunk number 11: eq2b
###################################################
out = fusedlasso1d(y)


###################################################
### code chunk number 12: eq2c
###################################################
plot(out, lambda=1)


###################################################
### code chunk number 13: fig2c
###################################################
plot(out, lambda=1)


###################################################
### code chunk number 14: eq8a
###################################################
set.seed(1)
y = matrix(runif(256), 16, 16)
i = (row(y) - 8.5)^2 + (col(y) - 8.5)^2 <= 4^2
y[i] = y[i] + 1


###################################################
### code chunk number 15: eq8b
###################################################
out = fusedlasso2d(y)


###################################################
### code chunk number 16: eq8x
###################################################
co = coef(out, nlam=5)


###################################################
### code chunk number 17: eq8c
###################################################
par(mar=c(1,1,2,1),mfrow=c(2,3))
cols = terrain.colors(30)
zlim = range(c(co$beta,y))
image(y,main=expression(y),col=cols,zlim=zlim,axes=FALSE)
for (i in 1:5) {
  image(matrix(co$beta[,i],nrow=16),col=cols,zlim=zlim,
  axes=FALSE)
  mtext(bquote(lambda==.(sprintf("%.3f",co$lambda[i]))))
}


###################################################
### code chunk number 18: fig8c
###################################################
par(mar=c(1,1,2,1),mfrow=c(2,3))
cols = terrain.colors(30)
zlim = range(c(co$beta,y))
image(y,main=expression(y),col=cols,zlim=zlim,axes=FALSE)
for (i in 1:5) {
  image(matrix(co$beta[,i],nrow=16),col=cols,zlim=zlim,
  axes=FALSE)
  mtext(bquote(lambda==.(sprintf("%.3f",co$lambda[i]))))
}


###################################################
### code chunk number 19: eq6a
###################################################
set.seed(1)
n = 100
i = 1:n
y = (i > 20 & i < 30) + 5*(i > 50 & i < 70) +
  rnorm(n, sd=0.1)
out = fusedlasso1d(y)
beta1 = coef(out, lambda=1.5)$beta


###################################################
### code chunk number 20: eq6b
###################################################
beta2 = softthresh(out, lambda=1.5, gamma=1)


###################################################
### code chunk number 21: eq7x
###################################################
plot(1:n, y, xlab="Position", ylab="Estimates")
abline(h=1.5, lty="dashed")
lines(1:n, beta1)
lines(1:n, beta2, col="red")
legend("topleft",lty=1,col=c("black","red"),
       legend=c(expression(gamma==0),expression(gamma==0.5)))


###################################################
### code chunk number 22: fig7x
###################################################
plot(1:n, y, xlab="Position", ylab="Estimates")
abline(h=1.5, lty="dashed")
lines(1:n, beta1)
lines(1:n, beta2, col="red")
legend("topleft",lty=1,col=c("black","red"),
       legend=c(expression(gamma==0),expression(gamma==0.5)))


###################################################
### code chunk number 23: eq3a
###################################################
set.seed(1)
n = 100
y = 50 - abs(1:n-n/2) + rnorm(n, sd=5)


###################################################
### code chunk number 24: eq3b
###################################################
out = trendfilter(y, ord=1)
plot(out, lambda=n)


###################################################
### code chunk number 25: fig3b
###################################################
out = trendfilter(y, ord=1)
plot(out, lambda=n)


###################################################
### code chunk number 26: eq4gib
###################################################
n = 100
y = sin(1:n/n*2*pi) + rnorm(n, sd=0.3)
out = trendfilter(y, ord=3)
plot(out, lambda=n)


###################################################
### code chunk number 27: fig4abc
###################################################
n = 100
y = sin(1:n/n*2*pi) + rnorm(n, sd=0.3)
out = trendfilter(y, ord=3)
plot(out, lambda=n)


###################################################
### code chunk number 28: eq5a
###################################################
set.seed(1)
n = 100
y = rep(sample(1:8,5), each=n/5) + rnorm(n, sd=0.8)
out = trendfilter(y, ord=0)


###################################################
### code chunk number 29: eq5b
###################################################
cv = cv.trendfilter(out)


###################################################
### code chunk number 30: eq5c
###################################################
plot(out, lambda=cv$lambda.min, main="Minimal CV error")


###################################################
### code chunk number 31: fig5c
###################################################
plot(out, lambda=cv$lambda.min, main="Minimal CV error")


###################################################
### code chunk number 32: eq5d
###################################################
plot(out, lambda=cv$lambda.1se, main="One standard error rule")


###################################################
### code chunk number 33: fig5d
###################################################
plot(out, lambda=cv$lambda.1se, main="One standard error rule")


