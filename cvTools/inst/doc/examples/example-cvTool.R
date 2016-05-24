library("robustbase")
data("coleman")
set.seed(1234)  # set seed for reproducibility

# set up function call for an MM regression model
call <- call("lmrob", formula = Y ~ .)
# set up folds for cross-validation
folds <- cvFolds(nrow(coleman), K = 5, R = 10)

# perform cross-validation
cvTool(call, data = coleman, y = coleman$Y, cost = rtmspe, 
    folds = folds, costArgs = list(trim = 0.1))
