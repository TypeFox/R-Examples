#### Test special cases  for  covMcd()

library(robustbase)

### 1) p = 1 ----------------------------------------------------
set.seed(1)
x <- c(rnorm(50),100, 1e10)
(r1 <- covMcd(x))
str(r1)
summary(r1)
## with alpha = 1
(r1.1 <- covMcd(x, alpha = 1))
str(r1.1)
summary(r1.1)

### 1b) p = 1, constant scale
(rc <- covMcd(rep(1,12)))
str(rc)
summary(rc)
## with alpha = 1
(rc1 <- covMcd(rep(1,12), alpha = 1))
str(rc1)
summary(rc1)

### 2)  constant observations  { multivariate scale == 0 } -----------
(X <- matrix(rep(2*(1:4), 12), nrow = 12, byrow = TRUE))
(rC  <- covMcd(X))
summary(rC)
(rC1 <- covMcd(X, alpha = 1))
summary(rC1)

### 3)  alpha = 1 : classical estimates --- for general cases --------


cat('Time elapsed: ', proc.time(),'\n') # for ``statistical reasons''
