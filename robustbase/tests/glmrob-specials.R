library(robustbase)

## Model without coefficients [ print.glmrob() tests for this ..]

### very simple model [with outliers]
set.seed(1)
y <- rpois(1000, lambda = 4)

## without outliers
m0o <- glm(y ~ 0, family = poisson, epsilon = 1e-12)
m1o <- glm(y ~ 1, family = poisson, epsilon = 1e-12)

y[1:3] <- 99:101 # outliers

m0 <- glm(y ~ 0, family = poisson, epsilon = 1e-12)
m1 <- glm(y ~ 1, family = poisson, epsilon = 1e-12)

## these both failed in version 0.1-2:
rm0 <- glmrob(y ~ 0, family = poisson, acc = 1e-12)
rm1 <- glmrob(y ~ 1, family = poisson, acc = 1e-12)

rm0
rm1
(s0 <- summary(rm0))
(s1 <- summary(rm1))
str(s1)
stopifnot(all.equal(c(coef(s1)),
		    c(1.390672035557, 0.016213613600955,
		      85.77187478275, 0), tolerance = 1e-13))# 32-b: 4.7e-15

cat('Time elapsed: ', proc.time(),'\n') # for ``statistical reasons''
