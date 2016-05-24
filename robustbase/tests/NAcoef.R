## test handing of NA coefficients / singular fits
## also check:
## -- what would have to be done if class "lm" was added.
## -- general compatibility to class lm.
require(robustbase)
source(system.file("test-tools-1.R", package="Matrix", mustWork=TRUE))
##-> assertError(), etc
ALL.equal <- function(x,y, ...) all.equal(x,y, tolerance=1e-15, ...)

## generate simple example data
data <- expand.grid(x1=letters[1:3], x2=LETTERS[1:3], rep=1:3)
set.seed(1)
data$y <- rnorm(nrow(data))
## drop all combinations of one interaction:
data <- subset(data, x1 != 'c' | (x2 != 'B' & x2 != 'C'))
## add collinear variables
data$x3 <- rnorm(nrow(data))
data$x4 <- rnorm(nrow(data))
data$x5 <- data$x3 + data$x4
## add some NA terms
data$y[1] <- NA
data$x4[2:3] <- NA ## to test anova

## Classical models start with 'cm', robust just with  'rm' (or just 'm'):
cm0 <- lm   (y ~ x1*x2 + x3,	       data)
cm1 <- lm   (y ~ x1*x2 + x3 + x4 + x5, data)
set.seed(2)
rm1 <- lmrob(y ~ x1*x2 + x3 + x4 + x5, data)
m3  <- lmrob(y ~ x1*x2 + x3 + x4,      data) # same column space as rm1
rm0 <- lmrob(y ~ x1*x2 + x3,	       data)

## clean version of rm1 (to check predict)
data2 <- data.frame(y=data$y[-(1:3)], rm1$x[,!is.na(rm1$coef)])
set.seed(2)
rm1c <- lmrob(y ~ x1b + x1c + x2B + x2C + x3 + x4 + x1b:x2B + x1b:x2C, data2)

## add class lm to rm1 (for now)
class(rm1) <- c(class(rm1), "lm")
class(rm0) <- c(class(rm0), "lm")

## the full matrix (data) should be returned by model matrix (frame)
stopifnot(all.equal(model.matrix(cm1), model.matrix(rm1)),
          all.equal(model.frame(cm1), model.frame(rm1)))
## qr decomposition should be for the full data and pivots identical lm result
qr.cm1 <- qr(cm1)$qr
qr.rm1 <- rm1$qr$qr
stopifnot(NCOL(qr.rm1) == NCOL(qr.cm1),
          NROW(qr.rm1) == NROW(qr.cm1),
          length(rm1$qr$qraux) == length(qr(cm1)$qraux),
          all.equal(rm1$qr$pivot, qr(cm1)$pivot),
          all.equal(dimnames(qr.rm1),dimnames(qr.cm1)))
## the alias function should return the same result
stopifnot(all.equal(alias(cm1), alias(rm1)))

####
## these helper functions should print NAs for the dropped coefficients
  print(rm1)
summary(rm1) -> s1
confint(rm1) -> ci1
stopifnot(identical(is.na(coef(cm1)), apply(ci1, 1L, anyNA)),
	  identical(sigma(rm1), s1$ sigma),
	  identical(vcov(rm1),  s1$ cov  ),
	  TRUE)

print(s1, showAlgo=FALSE)
ci1
## drop1 should return df = 0
#drop1(rm1) ## drop.lm does not return valid results (yet)!

####
## methods that should just drop the NA coefficients
## m3 is actually the same as rm1, so anova should raise an error
assertError(anova(rm1, m3, test="Wald"))
assertError(anova(rm1, m3, test="Deviance"))
## but comparing rm1 and rm0 should be ok
anova(rm1, rm0, test="Wald")
anova(rm1, rm0, test="Deviance")
## commands with single #:
## they do (or might) not return sensible results for robust fits
## and need to be checked again
#cooks.distance(rm1)
#deviance(rm1)
#dfbeta(rm1)
#dfbetas(rm1)
#effects(rm1) ## fails
#extractAIC(rm1)
#stopifnot(all.equal(hatvalues(rm1), robustbase:::lmrob.leverages(wqr=rm1$qr))) ## fails
#influence(rm1)
stopifnot(is.infinite(kr1 <- kappa(rm1)), kr1 == kappa(cm1), # = +Inf both
          identical(labels(rm1), labels(cm1)))

logLik(rm1)# well, and what does it mean?

## plot(rm1, which=1) ## plot.lmrob() fails "singular covariance" .. FIXME!
par(mfrow=c(2,2))
plot(rm1, which=2:4)
stopifnot(ALL.equal(predict(rm1), predict(rm1c)),
          ALL.equal(predict(rm1,  se.fit=TRUE, interval="confidence"),
		    predict(rm1c, se.fit=TRUE, interval="confidence")))
predict(rm1, type="terms", se.fit=TRUE, interval="confidence")
#proj(rm1) ## fails "FIXME"
residuals(rm1)
#rstandard(rm1)
#rstudent(rm1)
#simulate(rm1) ## just $weights needs to be changed to prior weights
V1 <- vcov(rm1) # but don't show the "eigen" part {vectors may flip sign}:
attributes(V1) <- attributes(V1)[c("dim","dimnames", "weights")]; V1
set.seed(12); sc <- simulate(cm1, 64)
set.seed(12); rc <- simulate(rm1, 64)

stopifnot(ALL.equal(sqrt(diag(V1)), coef(summary(rm1))[,"Std. Error"]),
	  all.equal(sc, rc, tolerance = 0.08),# dimension *and* approx. values (no NA)
	  identical(variable.names(rm1), variable.names(cm1)),
	  all.equal(residuals(rm1), residuals(cm1), tolerance = 0.05),# incl. names
	  all.equal(rstudent (rm1), rstudent (cm1), tolerance = 0.06),
	  identical(dimnames(rm1), dimnames(cm1)),
	  all.equal(dummy.coef(rm1), dummy.coef(cm1), tolerance= .5)) ## check mostly structure

## other helper functions
stopifnot(identical(case.names(rm1), case.names(cm1)),
          all.equal(family(rm1), family(cm1)),# identical() upto environment
          identical(formula(rm1), formula(cm1)),
          nobs(rm1) == nobs(cm1))
#add1(rm0, ~ . + x3 + x4 + x5) ## does not return valid results (yet)!


## test other initial estimators
lmrob(y ~ x1*x2 + x3 + x4 + x5, data, init="M-S")
lmrob(y ~ x1*x2 + x3 + x4 + x5, data, init=lmrob.lar)

## test all zero design matrix
data <- data.frame(y=1:10,x1=0,x2=0,os=2,w=c(0.5, 1))
(m5 <- lmrob(y ~ 1+x1+x2+offset(os), data, weights=w))
(sm5 <- summary(m5))
(m6 <- lmrob(y ~ 0+x1+x2+offset(os), data, weights=w))
(sm6 <- summary(m6))

sc5 <- summary(cm5 <- lm(y ~ 1+x1+x2+offset(os), data, weights=w))
sc6 <- summary(cm6 <- lm(y ~ 0+x1+x2+offset(os), data, weights=w))

stopifnot(all.equal(coef(m5), coef(cm5), tolerance = 0.01),
          identical(coef(m6), coef(cm6)),
          all.equal(coef(sm5), coef(sc5), tolerance = 0.05),
          identical(coef(sm6), coef(sc6)),
          identical(sm5$df, sc5$df),
          identical(sm6$df, sc6$df))
