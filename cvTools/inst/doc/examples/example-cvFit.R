library("robustbase")
data("coleman")

## via model fit
# fit an MM regression model
fit <- lmrob(Y ~ ., data=coleman)
# perform cross-validation
cvFit(fit, data = coleman, y = coleman$Y, cost = rtmspe, 
    K = 5, R = 10, costArgs = list(trim = 0.1), seed = 1234)

## via model fitting function
# perform cross-validation
# note that the response is extracted from 'data' in 
# this example and does not have to be supplied
cvFit(lmrob, formula = Y ~ ., data = coleman, cost = rtmspe, 
    K = 5, R = 10, costArgs = list(trim = 0.1), seed = 1234)

## via function call
# set up function call
call <- call("lmrob", formula = Y ~ .)
# perform cross-validation
cvFit(call, data = coleman, y = coleman$Y, cost = rtmspe, 
    K = 5, R = 10, costArgs = list(trim = 0.1), seed = 1234)
