data("coleman")
set.seed(1234)  # set seed for reproducibility

## set up folds for cross-validation
folds <- cvFolds(nrow(coleman), K = 5, R = 10)

## compare LS, MM and LTS regression

# perform cross-validation for an LS regression model
fitLm <- lm(Y ~ ., data = coleman)
cvLm <- perry(fitLm, splits = folds, 
    cost = rtmspe, trim = 0.1)

# perform cross-validation for an MM regression model
fitLmrob <- lmrob(Y ~ ., data = coleman, k.max = 500)
cvLmrob <- perry(fitLmrob, splits = folds, 
    cost = rtmspe, trim = 0.1)

# perform cross-validation for an LTS regression model
fitLts <- ltsReg(Y ~ ., data = coleman)
cvLts <- perry(fitLts, splits = folds, 
    cost = rtmspe, trim = 0.1)

# combine results into one object
cv <- perrySelect(LS = cvLm, MM = cvLmrob, LTS = cvLts)
cv

## convert MM regression results to data frame for plotting
# all replications for box plot
cvLmrobBox <- fortify(cvLmrob, reps = TRUE)
perryPlot(cvLmrobBox)
# aggregated results for dot plot
cvLmrobDot <- fortify(cvLmrob, reps = FALSE, seFactor = 1)
perryPlot(cvLmrobDot)

## convert combined results to data frame for plotting
# all replications for box plot
cvBox <- fortify(cv, reps = TRUE)
perryPlot(cvBox)
# aggregated results for dot plot
cvDot <- fortify(cv, reps = FALSE, seFactor = 1)
perryPlot(cvDot)
