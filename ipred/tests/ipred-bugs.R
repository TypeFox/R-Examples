library(ipred)


actversion <- paste(R.version$major, R.version$minor, sep=".")
thisversion <- "1.7.0"

#if (compareVersion(actversion, thisversion) >= 0) {
#  RNGversion("1.6.2")
#}
set.seed(29081975)

data("BreastCancer", package = "mlbench")
mod <- bagging(Class ~ Cl.thickness + Cell.size
                + Cell.shape + Marg.adhesion
                + Epith.c.size + Bare.nuclei
                + Bl.cromatin + Normal.nucleoli
                + Mitoses, data=BreastCancer, coob=TRUE)
print(mod)

print(a <- predict(mod, newdata=BreastCancer))
stopifnot(length(a) == nrow(BreastCancer))

# bagging failed if only one predictor was specified
# by Christoph M. Friedrich <chris@uni-wh.de>, April 29th, 2002

X <- as.data.frame(matrix(rnorm(1000), ncol=10))
y <- factor(ifelse(apply(X, 1, mean) > 0, 1, 0))
learn <- cbind(y, X)
mt <- bagging(y ~ V1, data=learn, coob=TRUE)
# <FIXME>
# This won't work because of some difficulties with predict.lda
# mt <- bagging(y ~ V1, data=learn, method="double", coob=FALSE)
# </FIXME>
X <- as.data.frame(matrix(rnorm(1000), ncol=10))
y <- apply(X, 1, mean) + rnorm(nrow(X))
learn <- cbind(y, X)
mt <- bagging(y ~ V1, data=learn, coob=TRUE)

# cv.numeric and bootest.numeric were broken, check for reasonaly values
X <- as.data.frame(matrix(rnorm(1000), ncol=10))
y <- apply(X, 1, mean) + rnorm(nrow(X))
learn <- cbind(y, X)
newy <- apply(X, 1, mean) + rnorm(nrow(X))
mod <- lm(y ~ ., data=learn)
trueerr <- sqrt(mean((newy - fitted(mod))^2))
cverr <- rep(0,5)
for (i in 1:5) cverr[i] <- errorest(y ~., data=learn, model=lm)$error
booterr <- errorest(y ~., data=learn, model=lm,
                    estimator="boot",est.para=control.errorest(nboot=50))$error
print(trueerr/mean(cverr))
print(trueerr/booterr)
