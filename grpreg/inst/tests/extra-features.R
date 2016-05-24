##source("~/dev/.grpreg.setup.R")
##require(grpreg)

##############################
.test = "logLik is correct" ##
##############################
n <- 50
group <- rep(0:4,5:1)
p <- length(group)
X <- matrix(rnorm(n*p),ncol=p)
y <- rnorm(n)
yy <- runif(n) > .5
fit.mle <- lm(y~X)
fit <- grpreg(X, y, group, penalty="grLasso", lambda.min=0)
check(logLik(fit)[100], logLik(fit.mle)[1], check.attributes=FALSE, tol=.001)
check(AIC(fit)[100], AIC(fit.mle), check.attributes=FALSE, tol=.001)
fit <- grpreg(X, y, group, penalty="gel", lambda.min=0)
check(logLik(fit)[100], logLik(fit.mle)[1], check.attributes=FALSE, tol=.001)
check(AIC(fit)[100], AIC(fit.mle), check.attributes=FALSE, tol=.001)
fit.mle <- glm(yy~X, family="binomial")
fit <- grpreg(X, yy, group, penalty="grLasso", lambda.min=0, family="binomial")
check(logLik(fit)[100], logLik(fit.mle)[1], check.attributes=FALSE, tol=.001)
check(AIC(fit)[100], AIC(fit.mle), check.attributes=FALSE, tol=.001)
fit <- grpreg(X, yy, group, penalty="gMCP", lambda.min=0, family="binomial")
check(logLik(fit)[100], logLik(fit.mle)[1], check.attributes=FALSE, tol=.001)
check(AIC(fit)[100], AIC(fit.mle), check.attributes=FALSE, tol=.001)
fit.mle <- glm(yy~X, family="poisson")
fit <- grpreg(X, yy, group, penalty="grLasso", lambda.min=0, family="poisson")
check(logLik(fit)[100], logLik(fit.mle)[1], check.attributes=FALSE, tol=.001)
check(AIC(fit)[100], AIC(fit.mle), check.attributes=FALSE, tol=.001)
fit <- grpreg(X, yy, group, penalty="gMCP", lambda.min=0, family="poisson")
check(logLik(fit)[100], logLik(fit.mle)[1], check.attributes=FALSE, tol=.001)
check(AIC(fit)[100], AIC(fit.mle), check.attributes=FALSE, tol=.001)

############################################
.test = "grpreg handles constant columns" ##
############################################
n <- 50
group <- rep(0:3,4:1)
p <- length(group)
X <- matrix(rnorm(n*p),ncol=p)
X[,group==2] <- 0
y <- rnorm(n)
yy <- y > 0
par(mfrow=c(3,3))
fit <- grpreg(X, y, group, penalty="grLasso"); plot(fit)
fit <- grpreg(X, y, group, penalty="gMCP"); plot(fit)
fit <- gBridge(X, y, group); plot(fit)
fit <- grpreg(X, yy, group, penalty="grLasso", family="binomial"); plot(fit)
fit <- grpreg(X, yy, group, penalty="gMCP", family="binomial"); plot(fit); fit$beta[,100]
fit <- gBridge(X, yy, group, family="binomial"); plot(fit); fit$beta[,100]
fit <- grpreg(X, yy, group, penalty="grLasso", family="poisson"); plot(fit)
fit <- grpreg(X, yy, group, penalty="gMCP", family="poisson"); plot(fit); fit$beta[,100]
fit <- gBridge(X, yy, group, family="poisson"); plot(fit); fit$beta[,100]

###################################################
.test = "grpreg handles groups of non-full rank" ##
###################################################
n <- 50
group <- rep(0:3,4:1)
p <- length(group)
X <- matrix(rnorm(n*p),ncol=p)
X[,7] <- X[,6]
y <- rnorm(n)
yy <- y > 0
par(mfrow=c(2,3))
fit <- grpreg(X, y, group, penalty="grLasso"); plot(fit)
fit <- grpreg(X, y, group, penalty="gMCP"); plot(fit)
fit <- grpreg(X, yy, group, penalty="grLasso", family="binomial"); plot(fit)
fit <- grpreg(X, yy, group, penalty="gMCP", family="binomial"); plot(fit)
fit <- grpreg(X, yy, group, penalty="grLasso", family="poisson"); plot(fit)
fit <- grpreg(X, yy, group, penalty="gMCP", family="poisson"); plot(fit)

#######################################
.test = "grpreg out-of-order groups" ##
#######################################
n <- 50
group <- rep(0:3,4:1)
ind <- sample(1:length(group))
p <- length(group)
X <- matrix(rnorm(n*p),ncol=p)
X[,group==2] <- 0
y <- rnorm(n)
yy <- y > 0
fit1 <- grpreg(X, y, group, penalty="grLasso")
fit2 <- grpreg(X[,ind], y, group[ind], penalty="grLasso")
b1 <- coef(fit1)[-1,][ind,]
b2 <- coef(fit2)[-1,]
check(b1, b2, tol=0.001, check.attributes=FALSE)

################################
.test = "grpreg named groups" ##
################################
n <- 50
group1 <- rep(0:3,4:1)
group2 <- rep(c("0", "A", "B", "C"), 4:1)
p <- length(group)
X <- matrix(rnorm(n*p),ncol=p)
X[,group==2] <- 0
y <- rnorm(n)
yy <- y > 0
fit1 <- grpreg(X, y, group1, penalty="grLasso")
fit2 <- grpreg(X, y, group2, penalty="grLasso")
check(coef(fit1), coef(fit2), tol=0.001, check.attributes=FALSE)

###################################################
.test = "cv.grpreg() seems to work" ##
###################################################
n <- 50
group <- rep(1:4,1:4)
p <- length(group)
X <- matrix(rnorm(n*p),ncol=p)
b <- rnorm(p, sd=2)
b[abs(b) < 1] <- 0
y <- rnorm(n, mean=X%*%b)
yy <- y > .5

par(mfrow=c(3,2))
require(glmnet)
cvfit <- cv.glmnet(X, y)
plot(cvfit)
cvfit <- cv.grpreg(X, y, group)
plot(cvfit)
cvfit <- cv.glmnet(X, yy, family="binomial")
plot(cvfit)
cvfit <- cv.grpreg(X, yy, group, family="binomial")
plot(cvfit)
cvfit <- cv.glmnet(X, yy, family="poisson")
plot(cvfit)
cvfit <- cv.grpreg(X, yy, group, family="poisson")
plot(cvfit)

###################################
.test = "group.multiplier works" ##
###################################
n <- 50
p <- 10
X <- matrix(rnorm(n*p),ncol=p)
y <- rnorm(n)
group <- rep(0:3,1:4)
gm <- 1:3
par(mfrow=c(3,2))
plot(fit <- grpreg(X, y, group, penalty="gMCP", lambda.min=0, group.multiplier=gm), main=fit$penalty)
plot(fit <- gBridge(X, y, group, lambda.min=0, group.multiplier=gm), main=fit$penalty)
plot(fit <- grpreg(X, y, group, penalty="grLasso", lambda.min=0, group.multiplier=gm), main=fit$penalty)
plot(fit <- grpreg(X, y, group, penalty="grMCP", lambda.min=0, group.multiplier=gm), main=fit$penalty)
plot(fit <- grpreg(X, y, group, penalty="grSCAD", lambda.min=0, group.multiplier=gm), main=fit$penalty)

##################################################
.test = "cv.grpreg() options work for gaussian" ##
##################################################
n <- 50
group <- rep(0:4,5:1)
p <- length(group)
X <- matrix(rnorm(n*p),ncol=p)
b <- c(-3, 3, rep(0, 13))
y <- rnorm(n, mean=X%*%b, sd=1)

par(mfrow=c(2,2))
cvfit <- cv.grpreg(X, y, group=group)
plot(cvfit, type="all")
summary(cvfit)

b <- c(-3, 3, rep(0, 13))
y <- rnorm(n, mean=X%*%b, sd=5)
cvfit <- cv.grpreg(X, y, group=group)
plot(cvfit, type="all")
coef(cvfit)
predict(cvfit, type="coef")
predict(cvfit, type="vars")
predict(cvfit, type="groups")
predict(cvfit, type="norm")
predict(cvfit, X)

b <- rep(0, 15)
y <- rnorm(n, mean=X%*%b, sd=5)
cvfit <- cv.grpreg(X, y, group=group)
plot(cvfit, type="all")

##################################################
.test = "cv.grpreg() options work for binomial" ##
##################################################
n <- 200
group <- rep(0:4,5:1)
p <- length(group)
X <- matrix(rnorm(n*p),ncol=p)
b <- c(-3, 3, rep(0, 13))
y <- rnorm(n, mean=X%*%b, sd=1) > 0.5

par(mfrow=c(2,2))
cvfit <- cv.grpreg(X, y, group=group, family="binomial")
plot(cvfit, type="all")
summary(cvfit)

b <- c(-3, 3, rep(0, 13))
y <- rnorm(n, mean=X%*%b, sd=5) > 0.5
cvfit <- cv.grpreg(X, y, group=group, family="binomial")
plot(cvfit, type="all")
coef(cvfit)
predict(cvfit, type="coef")
predict(cvfit, type="vars")
predict(cvfit, type="groups")
predict(cvfit, type="norm")
predict(cvfit, X)
predict(cvfit, X, type="response")

b <- rep(0, 15)
y <- rnorm(n, mean=X%*%b, sd=5) > 0.5
cvfit <- cv.grpreg(X, y, group=group, family="binomial")
plot(cvfit, type="all")

##################################################
.test = "cv.grpreg() options work for poisson" ##
##################################################
n <- 100
group <- rep(1:4, 4:1)
p <- length(group)
X <- matrix(rnorm(n*p),ncol=p)
b <- c(-3, 3, rep(0, p-2))
y <- floor(exp(X%*%b + rnorm(n)))

par(mfrow=c(3,1))
cvfit <- cv.grpreg(X, y, group=group, family="poisson")
plot(cvfit, type="all")
summary(cvfit)

coef(cvfit)
predict(cvfit, type="coef")
predict(cvfit, type="vars")
predict(cvfit, type="groups")
predict(cvfit, type="norm")
predict(cvfit, X)
predict(cvfit, X, type="response")

y <- sample(0:5, n, replace=TRUE)
cvfit <- cv.grpreg(X, y, group=group, family="poisson")
plot(cvfit, type="all")

#############################
.test = "dfmax, gmax work" ##
#############################
n <- 100
group <- rep(1:10, rep(3,10))
p <- length(group)
X <- matrix(rnorm(n*p),ncol=p)
y <- rnorm(n)
yy <- runif(n) > .5

## dfmax
dfmax <- 21
fit <- grpreg(X, y, group, penalty="grLasso", lambda.min=0, dfmax=dfmax)
nv <- sapply(predict(fit, type="vars"), length)
check(max(head(nv, length(nv)-1)) <= dfmax)
check(max(nv) > 3)
fit <- grpreg(X, y, group, penalty="gel", lambda.min=0, dfmax=dfmax)
nv <- sapply(predict(fit, type="vars"), length)
check(max(head(nv, length(nv)-1)) <= dfmax)
check(max(nv) > 3)
fit <- grpreg(X, yy, group, penalty="grLasso", family="binomial", lambda.min=0, dfmax=dfmax)
nv <- sapply(predict(fit, type="vars"), length)
check(max(head(nv, length(nv)-1)) <= dfmax)
check(max(nv) > 3)
fit <- grpreg(X, yy, group, penalty="gel", family="binomial", lambda.min=0, dfmax=dfmax)
nv <- sapply(predict(fit, type="vars"), length)
check(max(head(nv, length(nv)-1)) <= dfmax)
check(max(nv) > 3)

## gmax
gmax <- 7
fit <- grpreg(X, y, group, penalty="grLasso", lambda.min=0, gmax=gmax)
ng <- sapply(predict(fit, type="groups"), length)
check(max(head(ng, length(ng)-1)) <= gmax)
check(max(ng) > 2)
fit <- grpreg(X, y, group, penalty="gel", lambda.min=0, gmax=gmax)
ng <- sapply(predict(fit, type="groups"), length)
check(max(head(ng, length(ng)-1)) <= gmax)
check(max(ng) > 2)
fit <- grpreg(X, yy, group, penalty="grLasso", family="binomial", lambda.min=0, gmax=gmax)
ng <- sapply(predict(fit, type="groups"), length)
check(max(head(ng, length(ng)-1)) <= gmax)
check(max(ng) > 2)
fit <- grpreg(X, yy, group, penalty="gel", family="binomial", lambda.min=0, gmax=gmax)
ng <- sapply(predict(fit, type="groups"), length)
check(max(head(ng, length(ng)-1)) <= gmax)
check(max(ng) > 2)
