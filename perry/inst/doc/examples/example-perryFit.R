data("coleman")
set.seed(1234)  # set seed for reproducibility

## via model fit
# fit an MM regression model
fit <- lmrob(Y ~ ., data=coleman)
# perform cross-validation
perryFit(fit, data = coleman, y = coleman$Y, 
    splits = foldControl(K = 5, R = 10), 
    cost = rtmspe, costArgs = list(trim = 0.1), 
    seed = 1234)

## via model fitting function
# perform cross-validation
# note that the response is extracted from 'data' in 
# this example and does not have to be supplied
perryFit(lmrob, formula = Y ~ ., data = coleman, 
    splits = foldControl(K = 5, R = 10), 
    cost = rtmspe, costArgs = list(trim = 0.1), 
    seed = 1234)

## via function call
# set up function call
call <- call("lmrob", formula = Y ~ .)
# perform cross-validation
perryFit(call, data = coleman, y = coleman$Y, 
    splits = foldControl(K = 5, R = 10), 
    cost = rtmspe, costArgs = list(trim = 0.1), 
    seed = 1234)
