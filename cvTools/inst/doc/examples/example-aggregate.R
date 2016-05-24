library("robustbase")
data("coleman")
set.seed(1234)  # set seed for reproducibility

## set up folds for cross-validation
folds <- cvFolds(nrow(coleman), K = 5, R = 10)


## compare raw and reweighted LTS estimators for 
## 50% and 75% subsets

# 50% subsets
fitLts50 <- ltsReg(Y ~ ., data = coleman, alpha = 0.5)
cvFitLts50 <- cvLts(fitLts50, cost = rtmspe, folds = folds, 
    fit = "both", trim = 0.1)

# 75% subsets
fitLts75 <- ltsReg(Y ~ ., data = coleman, alpha = 0.75)
cvFitLts75 <- cvLts(fitLts75, cost = rtmspe, folds = folds, 
    fit = "both", trim = 0.1)

# combine results into one object
cvFitsLts <- cvSelect("0.5" = cvFitLts50, "0.75" = cvFitLts75)
cvFitsLts

# summary of the results with the 50% subsets
aggregate(cvFitLts50, summary)
# summary of the combined results
aggregate(cvFitsLts, summary)


## evaluate MM regression models tuned for 
## 80%, 85%, 90% and 95% efficiency
tuning <- list(tuning.psi=c(3.14, 3.44, 3.88, 4.68))

# set up function call
call <- call("lmrob", formula = Y ~ .)
# perform cross-validation
cvFitsLmrob <- cvTuning(call, data = coleman, 
    y = coleman$Y, tuning = tuning, cost = rtmspe, 
    folds = folds, costArgs = list(trim = 0.1))
cvFitsLmrob

# summary of results
aggregate(cvFitsLmrob, summary)
