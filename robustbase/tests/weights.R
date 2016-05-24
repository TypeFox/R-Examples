## test handing of weights and offset argument
require(robustbase)

## generate simple example data
data <- expand.grid(x1=letters[1:3], x2=LETTERS[1:4], rep=1:3)
## generate offset column
data$os <- 1:nrow(data)
set.seed(1)
data$y <- data$os + rnorm(nrow(data))
## add collinear variables
data$x3 <- rnorm(nrow(data))
data$x4 <- rnorm(nrow(data))
data$x5 <- data$x3 + data$x4
## add some NA terms
data$y[1] <- NA
data$x4[2:3] <- NA ## to test anova
## generate weights
## some obs with weight 0
data$weights <- as.numeric(with(data, x1 != 'c' | (x2 != 'B' & x2 != 'C')))
## some obs with weight 2
data$weights[data$x1 == 'b'] <- 2
data2 <- rbind(subset(data, weights>0), subset(data, weights==2))

## using these parameters we're essentially forcing lmrob() to
## fit a classic model --> easier to compare to lm()
ctrl <- lmrob.control(psi="optimal", tuning.chi = 20, bb = 0.0003846154,
                      tuning.psi=20, method="SM", cov=".vcov.w")

## Classical models start with 'cm', robust just with  'rm' (or just 'm'):
(cm0 <- lm   (y ~ x1*x2 + x3 + x4 + x5 + offset(os), data))
(cm1 <- lm   (y ~ x1*x2 + x3 + x4 + x5 + offset(os), data,  weights=weights))
(cm2 <- lm   (y ~ x1*x2 + x3 + x4 + x5,              data2, offset=os))
(rm0 <- lmrob(y ~ x1*x2 + x3 + x4 + x5 + offset(os), data,                   control=ctrl))
set.seed(2)
(rm1 <- lmrob(y ~ x1*x2 + x3 + x4 + x5 + offset(os), data,  weights=weights, control=ctrl))
set.seed(2)
(rm2 <- lmrob(y ~ x1*x2 + x3 + x4 + x5,              data2, offset=os,       control=ctrl))

sc0 <- summary(cm0)
sc1 <- summary(cm1)
sc2 <- summary(cm2)
(sr0 <- summary(rm0))
(sr1 <- summary(rm1))
(sr2 <- summary(rm2))

## test Estimates, Std. Errors, ...
stopifnot(all.equal(coef(cm1), coef(cm2)),
          all.equal(coef(rm1), coef(rm2)),
          all.equal(coef(sc0), coef(sr0)),
          all.equal(coef(sc1), coef(sr1)),
          all.equal(coef(sc2), coef(sr2)))

## test class "lm" methods that do not depend on weights
meths1 <- c("family",
            "formula",
            "labels",
            "model.matrix",
            "na.action",
            "terms")
for (meth in meths1)
    stopifnot(all.equal(do.call(meth, list(rm0)),
                        do.call(meth, list(rm1))))

## class "lm" methods that depend on weights
##                                      FIXME:
meths2 <- c(#"AIC",
            "alias",
            #"BIC",
            "case.names",
            "coef",
            "confint",
            #"cooks.distance",
            #"deviance",
            "df.residual",
            #"dfbeta",
            #"dfbetas",
            #"drop1",
            "dummy.coef",
            #"effects",
            #"extractAIC",
            #"hatvalues",
            #"influence",
            "kappa",
            #"logLik",
            #"model.frame", ## disable because of zero.weights attribute
            "nobs",
            "predict",
                                        #"proj",
                                        #"rstandard",
                                        #"rstudent",
                                        #"simulate",
            ##"summary", ## see above
            "variable.names",
            ##"vcov",    ## see below
            "weights")
op <- options(warn = 1)# print immediately
for (meth in meths2) {
    cat(meth,":")
    .SW. <- if(meth == "weights") suppressWarnings else identity # for suppressing
    ## No weights defined for this object. Use type="robustness" ....
    stopifnot(all.equal(do.call(meth, list(cm1)),
                        do.call(meth, list(rm1))),
              all.equal(do.call(meth, list(cm2)),
		   .SW.(do.call(meth, list(rm2)))))

    cat("\n")
}
options(op)# reverting

## further tests:
anova(rm1, update(rm1, ~ . - x4 - x5))
anova(rm2, update(rm2, ~ . - x4 - x5))

stopifnot(all.equal(fitted(cm0),          fitted(rm0)),
          all.equal(fitted(cm1),          fitted(rm1)),
          ## FIXME?: fitted(cm2) is of class AsIs but fitted(rm2) is numeric
          all.equal(unclass(fitted(cm2)), fitted(rm2)))

nd <- expand.grid(x1=letters[1:3], x2=LETTERS[1:4])
set.seed(3)
nd$x3 <- rnorm(nrow(nd))
nd$x4 <- rnorm(nrow(nd))
nd$x5 <- rnorm(nrow(nd))
nd$os <- nrow(nd):1
wts   <- runif(nrow(nd))
stopifnot(all.equal(predict(cm0, nd, interval="prediction"),
                    predict(rm0, nd, interval="prediction")),
          all.equal(predict(cm1, nd, interval="prediction"),
                    predict(rm1, nd, interval="prediction")),
          all.equal(predict(cm2, nd, interval="prediction"),
                    predict(rm2, nd, interval="prediction")),
          all.equal(predict(cm0, nd, interval="prediction", weights=wts),
                    predict(rm0, nd, interval="prediction", weights=wts)),
          all.equal(predict(cm1, nd, interval="prediction", weights=wts),
                    predict(rm1, nd, interval="prediction", weights=wts)),
          all.equal(predict(cm2, nd, interval="prediction", weights=wts),
                    predict(rm2, nd, interval="prediction", weights=wts),
                    tolerance=1e-7))

## Padding can lead to differing values here
## so test only full rank part
qrEQ <- function(m1, m2) {
    q1 <- qr(m1)
    q2 <- qr(m2)
    r <- 1:q1$rank
    stopifnot(q1$rank == q2$rank,
              all.equal(q1$pivot, q2$pivot),
              all.equal(q1$qraux[r],q2$qraux[r]),
              all.equal(q1$qr[r,r], q2$qr[r,r]))
}
qrEQ(cm0, rm0)
qrEQ(cm1, rm1)
qrEQ(cm2, rm2)

stopifnot(all.equal(residuals(cm0),                      residuals(rm0)),
          all.equal(residuals(cm1),                      residuals(rm1)),
          ## FIXME?: residuals(cm2) is of class AsIs but residuals(rm2) is numeric
          all.equal(unclass(residuals(cm2)),             residuals(rm2)),
          all.equal(resid(cm0, type="pearson"),          resid(rm0, type="pearson")),
          all.equal(resid(cm1, type="pearson"),          resid(rm1, type="pearson")),
          all.equal(unclass(resid(cm2, type="pearson")), resid(rm2, type="pearson")))

stopifnot(all.equal(vcov(cm0), vcov(rm0), check.attributes=FALSE),
          all.equal(vcov(cm1), vcov(rm1), check.attributes=FALSE),
          all.equal(vcov(cm2), vcov(rm2), check.attributes=FALSE))

## Null fits (rank(X)==0) are tested in NAcoef.R

## testing weight=0 bug
lmrob(y ~ x3, data, weights=weights)
