data("coleman")

# perform cross-validation for an LTS regression model
fit <- ltsReg(Y ~ ., data = coleman)
folds <- foldControl(K = 5, R = 10)
cv <- perry(fit, splits = folds, fit = "both", 
    cost = rtmspe, trim = 0.1, seed = 1234)

# compare original and reshaped object
cv
perryReshape(cv)
