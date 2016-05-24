## Tests concerning seemingly unrelated regressions/multitask learning
##source("~/dev/.grpreg.setup.R")

#####################################
.test = "multitask learning works" ##
#####################################
n <- 50
p <- 10
X <- matrix(rnorm(n*p),ncol=p)
b <- c(-3, 3, -1, 1, rep(0, 6))
Y <- cbind(rnorm(n, mean=X%*%b, sd=1),
           5+rnorm(n, mean=X%*%b, sd=1),
           10+rnorm(n, mean=X%*%b, sd=1))
colnames(X) <- LETTERS[1:10]

par(mfcol=c(3,2))
fit <- grpreg(X, Y, penalty="grLasso"); plot(fit)
fit <- grpreg(X, Y, penalty="cMCP"); plot(fit)
fit <- gBridge(X, Y); plot(fit)

n <- 200
X <- matrix(rnorm(n*p),ncol=p)
Y <- cbind(rnorm(n, mean=X%*%b/10, sd=1),
           rnorm(n, mean=X%*%b/10, sd=1),
           rnorm(n, mean=X%*%b/10, sd=1)) > 0
fit <- grpreg(X, Y, penalty="grLasso", family="binomial", lambda.min=0.1); plot(fit)
fit <- grpreg(X, Y, penalty="cMCP", family="binomial", lambda.min=0.2); plot(fit)
fit <- gBridge(X, Y, family="binomial", lambda.min=0.2); plot(fit)

#####################################################
.test = "coef/predict work for multitask learning" ##
#####################################################
n <- 50
p <- 10
X <- matrix(rnorm(n*p),ncol=p)
b <- c(-3, 3, -1, 1, rep(0, 6))
Y <- cbind(rnorm(n, mean=X%*%b, sd=1),
           5+rnorm(n, mean=X%*%b, sd=1),
           10+rnorm(n, mean=X%*%b, sd=1))
colnames(X) <- LETTERS[1:10]

fit <- grpreg(X, Y)
coef(fit, which=1:2)
coef(fit, lambda=1)
predict(fit, lambda=1, type="nvars")
predict(fit, which=c(30,60), type="nvars")
predict(fit, lambda=1, type="ngroups")
predict(fit, which=c(30,60), type="ngroups")
predict(fit, lambda=1, type="groups")
predict(fit, which=c(30,60), type="groups")
predict(fit, lambda=1, type="norm")
predict(fit, which=c(30,60), type="norm")
head(predict(fit, X, lambda=1))
predict(fit, X, which=c(30,60))[1:10,,]

n <- 200
X <- matrix(rnorm(n*p),ncol=p)
Y <- cbind(rnorm(n, mean=X%*%b/3, sd=1),
           rnorm(n, mean=X%*%b/3, sd=1),
           rnorm(n, mean=X%*%b/3, sd=1)) > 0
fit <- grpreg(X, Y, family="binomial")
predict(fit, lambda=0.1, type="nvars")
predict(fit, lambda=0.1, type="ngroups")
predict(fit, lambda=0.1, type="groups")
predict(fit, lambda=0.1, type="norm")
head(predict(fit, X, lambda=0.1))
head(predict(fit, X, lambda=0.1, type="response"))
head(predict(fit, X, lambda=0.1, type="class"))

############################################################
.test = "multitask learning reproduces linear regression" ##
############################################################
n <- 50
p <- 10
X <- matrix(rnorm(n*p),ncol=p)
Y <- matrix(rnorm(n*3),ncol=3)
fit.mle <- lm(Y~X)
reg <- coef(fit.mle)
cMCP <- coef(fit <- grpreg(X, Y, penalty="cMCP", lambda.min=0), which=100)
check(t(cMCP), reg, tolerance=0.01, check.attributes=FALSE)
p <- predict(fit, X, which=100)
check(p, predict(fit.mle), tolerance=.01, check.attributes=FALSE)

bridge <- coef(fit <- gBridge(X, Y, lambda.min=0), which=1)
check(t(bridge), reg, tolerance=.01, check.attributes=FALSE)
p <- predict(fit, X, which=1)
check(p, predict(fit.mle), tolerance=.01, check.attributes=FALSE)

grLasso <- coef(fit <- grpreg(X, Y, penalty="grLasso", lambda.min=0), which=100)
check(t(grLasso), reg, tolerance=.01, check.attributes=FALSE)
p <- predict(fit, X, which=100)
check(p, predict(fit.mle), tolerance=.01, check.attributes=FALSE)

grMCP <- coef(fit <- grpreg(X, Y, penalty="grMCP", lambda.min=0), which=100)
check(t(grMCP), reg, tolerance=.01, check.attributes=FALSE)
p <- predict(fit, X, which=100)
check(p, predict(fit.mle), tolerance=.01, check.attributes=FALSE)

grSCAD <- coef(fit <- grpreg(X, Y, penalty="grSCAD", lambda.min=0), which=100)
check(t(grSCAD), reg, tolerance=.01, check.attributes=FALSE)
p <- predict(fit, X, which=100)
check(p, predict(fit.mle), tolerance=.01, check.attributes=FALSE)

##############################################################
.test = "multitask learning reproduces logistic regression" ##
##############################################################
n <- 200
p <- 10
X <- matrix(rnorm(n*p),ncol=p)
Y <- matrix(rnorm(n*3),ncol=3)>0
fit.mle <- glm(Y[,3]~X, family=binomial)
mle <- coef(fit.mle)

beta <- coef(fit <- grpreg(X, Y, lambda.min=0, family="binomial"), which=100)[3,]
check(beta, mle, tolerance=.01, check.attributes=FALSE)
p <- predict(fit, X, which=100)[,3]
check(p, predict(fit.mle), tolerance=.01, check.attributes=FALSE)
p <- predict(fit, X, which=100, type="response")[,3]
check(p, predict(fit.mle, type="response"), tolerance=.01, check.attributes=FALSE)

bridge <- coef(fit <- gBridge(X, Y, family="binomial", lambda.min=0), which=1)[3,]
check(bridge, mle, tolerance=.01, check.attributes=FALSE)
p <- predict(fit, X, which=1)[,3]
check(p, predict(fit.mle), tolerance=.01, check.attributes=FALSE)
p <- predict(fit, X, which=1, type="response")[,3]
check(p, predict(fit.mle, type="response"), tolerance=.01, check.attributes=FALSE)

##########################################################
.test = "cross-validation for multitask learning works" ##
##########################################################
n <- 50
p <- 10
X <- matrix(rnorm(n*p),ncol=p)
colnames(X) <- LETTERS[1:10]
b <- c(-3, 3, -1, 1, rep(0, 6))
Y <- cbind(rnorm(n, mean=X%*%b, sd=1),
           5+rnorm(n, mean=X%*%b, sd=1),
           10+rnorm(n, mean=X%*%b, sd=1))

par(mfcol=c(2,2))
cvfit <- cv.grpreg(X, Y); plot(cvfit, type="all"); summary(cvfit)
cvfit <- cv.grpreg(X, Y, penalty="cMCP"); plot(cvfit, type="all"); summary(cvfit)

b <- rep(0,10)
Y <- cbind(rnorm(n, mean=X%*%b, sd=1),
           5+rnorm(n, mean=X%*%b, sd=1),
           10+rnorm(n, mean=X%*%b, sd=1))
cvfit <- cv.grpreg(X, Y); plot(cvfit, type="all"); summary(cvfit)
cvfit <- cv.grpreg(X, Y, penalty="cMCP"); plot(cvfit, type="all"); summary(cvfit)

n <- 200
X <- matrix(rnorm(n*p),ncol=p)
b <- c(-3, 3, -1, 1, rep(0, 6))
Y <- cbind(rnorm(n, mean=X%*%b, sd=1),
           rnorm(n, mean=X%*%b, sd=1),
           rnorm(n, mean=X%*%b, sd=1)) > 0
cvfit <- cv.grpreg(X, Y, family="binomial"); plot(cvfit, type="all"); summary(cvfit)
cvfit <- cv.grpreg(X, Y, penalty="cMCP", family="binomial"); plot(cvfit, type="all"); summary(cvfit)

b <- rep(0,10)
Y <- cbind(rnorm(n, mean=X%*%b, sd=1),
           rnorm(n, mean=X%*%b, sd=1),
           rnorm(n, mean=X%*%b, sd=1)) > 0
cvfit <- cv.grpreg(X, Y, family="binomial"); plot(cvfit, type="all"); summary(cvfit)
cvfit <- cv.grpreg(X, Y, penalty="cMCP", family="binomial"); plot(cvfit, type="all"); summary(cvfit)
