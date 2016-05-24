data("coleman")
set.seed(1234)  # set seed for reproducibility

## set up folds for cross-validation
folds <- cvFolds(nrow(coleman), K = 5, R = 10)

## compare raw and reweighted LTS estimators for 50% and 75% 
## subsets based on their RTMSPE with 25% trimming

# 50% subsets
fit50 <- ltsReg(Y ~ ., data = coleman, alpha = 0.5)
cv50 <- perry(fit50, splits = folds, fit = "both", 
    cost = rtmspe, trim = 0.25)

# 75% subsets
fit75 <- ltsReg(Y ~ ., data = coleman, alpha = 0.75)
cv75 <- perry(fit75, splits = folds, fit = "both", 
    cost = rtmspe, trim = 0.25)

# combine results into one object
cv <- perrySelect("0.5" = cv50, "0.75" = cv75)
cv

## recompute the RTMSPE with 10% trimming
reperry(cv50, cost = rtmspe, trim = 0.1)
reperry(cv, cost = rtmspe, trim = 0.1)
