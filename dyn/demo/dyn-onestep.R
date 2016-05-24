
# one step ahead prediction

# test data
set.seed(1)
y <- ts(rnorm(10, 1:10, 0.1))

# fit model
library(dyn)
y.lm <- dyn$lm(y ~ lag(y,-1))

# use predict
tail(predict(y.lm, y), 1)

# or multiply by coefficients giving same result
coef(y.lm) %*% c(1, tail(y,1))

# Now try it using quantile regression
library(quantreg)
y.rq <- dyn$rq(y ~ lag(y,-1))
tail(predict(y.rq, y), 1)
coef(y.rq) %*% c(1, tail(y,1))

