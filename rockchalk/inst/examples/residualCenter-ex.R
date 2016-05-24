set.seed(123)
x1 <- rnorm(100)
x2 <- rnorm(100)
x3 <- rnorm(100)
x4 <- rnorm(100)
y <- rnorm(100)
dat <- data.frame(y, x1,x2,x3,x4)
rm(x1,x2,x3,x4,y)
m1 <- lm(y~ x1*x2 + x4, data = dat)
 
m1RC <- residualCenter(m1)

m1RCs <- summary(m1RC)
## The stage 1 centering regressions can be viewed as well
## lapply(m1RC$rcRegressions, summary)

## Verify residualCenter() output against the manual calculation
dat$x1rcx2 <- as.numeric(resid(lm(I(x1*x2) ~ x1 + x2, data = dat)))
m1m <- lm(y ~ x1 + x2 + x4 + x1rcx2, data=dat)
summary(m1m)
cbind("residualCenter" = coef(m1RC), "manual" = coef(m1m))


m2 <- lm(y~ x1*x2*x3 + x4, data=dat)
m2RC <- residualCenter(m2)
m2RCs <- summary(m2RC)

## Verify that result manually
dat$x2rcx3 <- as.numeric(resid(lm(I(x2*x3) ~ x2 + x3, data = dat)))
dat$x1rcx3 <- as.numeric(resid(lm(I(x1*x3) ~ x1 + x3, data = dat)))
dat$x1rcx2rcx3 <- as.numeric( resid(lm(I(x1*x2*x3) ~ x1 + x2 + x3 + x1rcx2 +
                                       x1rcx3 + x2rcx3 , data=dat)))
(m2m <- lm(y ~ x1 + x2 + x3+ x4 + x1rcx2 + x1rcx3 + x2rcx3 + x1rcx2rcx3,
           data = dat))

cbind("residualCenter" = coef(m2RC), "manual" = coef(m2m))


### As good as pequod's lmres
### not run because pequod generates R warnings
###
### if (require(pequod)){
###  pequodm1 <- lmres(y ~ x1*x2*x3 + x4, data=dat) 
###  pequodm1s <- summary(pequodm1)
###  coef(pequodm1s)
### }

### Works with any number of interactions. See:

m3 <- lm(y~ x1*x2*x3*x4, data=dat)
m3RC <- residualCenter(m3)
summary(m3RC)
##'
## Verify that one manually (Gosh, this is horrible to write out)
dat$x1rcx4 <- as.numeric(resid(lm(I(x1*x4) ~ x1 + x4, data=dat)))
dat$x2rcx4 <- as.numeric(resid(lm(I(x2*x4) ~ x2 + x4, data=dat)))
dat$x3rcx4 <- as.numeric(resid(lm(I(x3*x4) ~ x3 + x4, data=dat)))
dat$x1rcx2rcx4 <- as.numeric(resid(lm(I(x1*x2*x4) ~ x1 + x2 + x4 +
                                      x1rcx2 + x1rcx4 + x2rcx4, data=dat)))
dat$x1rcx3rcx4 <- as.numeric(resid(lm(I(x1*x3*x4) ~ x1 + x3 + x4 +
                                      x1rcx3 + x1rcx4 + x3rcx4, data=dat)))
dat$x2rcx3rcx4 <- as.numeric(resid(lm(I(x2*x3*x4) ~ x2 + x3 + x4 +
                                      x2rcx3 + x2rcx4 + x3rcx4, data=dat)))
dat$x1rcx2rcx3rcx4 <-
    as.numeric(resid(lm(I(x1*x2*x3*x4) ~ x1 + x2 + x3 + x4 +
                        x1rcx2 + x1rcx3 + x2rcx3 + x1rcx4  + x2rcx4 +
                        x3rcx4  + x1rcx2rcx3 + x1rcx2rcx4 + x1rcx3rcx4 +
                        x2rcx3rcx4, data=dat)))
(m3m <- lm(y ~ x1 + x2 + x3 + x4 + x1rcx2 + x1rcx3 + x2rcx3 + x1rcx4 +
           x2rcx4 + x3rcx4 + x1rcx2rcx3 + x1rcx2rcx4 + x1rcx3rcx4 +
           x2rcx3rcx4 + x1rcx2rcx3rcx4, data=dat))

cbind("residualCenter"=coef(m3RC), "manual"=coef(m3m))

### If you want to fit a sequence of models, as in pequod, can do.

tm <-terms(m2)
tmvec <- attr(terms(m2), "term.labels")
f1 <- tmvec[grep(":", tmvec, invert = TRUE)]
f2 <- tmvec[grep(":.*:", tmvec, invert = TRUE)]
f3 <- tmvec[grep(":.*:.*:", tmvec, invert = TRUE)]

## > f1
## [1] "x1" "x2" "x3" "x4"
## > f2
## [1] "x1"    "x2"    "x3"    "x4"    "x1:x2" "x1:x3" "x2:x3"
## > f3
## [1] "x1"       "x2"       "x3"       "x4"       "x1:x2"    "x1:x3"    "x2:x3"   
## [8] "x1:x2:x3"

f1 <- lm(as.formula(paste("y","~", paste(f1, collapse=" + "))), data=dat)
f1RC <- residualCenter(f1)
summary(f1RC)

f2 <- lm(as.formula(paste("y","~", paste(f2, collapse=" + "))), data=dat)
f2RC <- residualCenter(f2)
summary(f2RC)

f3 <- lm(as.formula(paste("y","~", paste(f3, collapse=" + "))), data=dat)
f3RC <- residualCenter(f3)
summary(f3RC)
