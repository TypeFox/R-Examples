library("robustbase")
data("coleman")

## evaluate MM regression models tuned for 85% and 95% efficiency
tuning <- list(tuning.psi = c(3.443689, 4.685061))

## via model fitting function
# perform cross-validation
# note that the response is extracted from 'data' in 
# this example and does not have to be supplied
cvTuning(lmrob, formula = Y ~ ., data = coleman, tuning = tuning, 
    cost = rtmspe, K = 5, R = 10, costArgs = list(trim = 0.1), 
    seed = 1234)

## via function call
# set up function call
call <- call("lmrob", formula = Y ~ .)
# perform cross-validation
cvTuning(call, data = coleman, y = coleman$Y, tuning = tuning, 
    cost = rtmspe, K = 5, R = 10, costArgs = list(trim = 0.1), 
    seed = 1234)
