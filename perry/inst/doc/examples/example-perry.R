## load data and fit an LS regression model
data("mtcars")
fit <- lm(mpg ~ wt + cyl, data=mtcars)

## perform cross-validation
# K-fold CV
perry(fit, foldControl(K = 5, R = 10), seed = 1234)
# leave-one-out CV
perry(fit, foldControl(K = nrow(mtcars)))

## perform random splitting
perry(fit, splitControl(m = 6, R = 10), seed = 1234)

## perform bootstrap prediction error estimation
# 0.632 estimator
perry(fit, bootControl(R = 10, type = "0.632"), seed = 1234)
# out-of-bag estimator
perry(fit, bootControl(R = 10, type = "out-of-bag"), seed = 1234)
