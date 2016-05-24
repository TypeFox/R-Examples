set.seed(20121204)

# generate some balance data
nvars <- 10
varnames <- paste("V", letters[1:nvars])

balance_data <- matrix(c(rnorm(n = nvars, mean = 1, sd = 0.5), 
                         rnorm(n = nvars, mean = 0, sd = 0.5)),
                       ncol = 2)

colnames(balance_data) <- c("Before Adjustment", "After Matching")

rownames(balance_data) <- varnames

RItools:::balanceplot(balance_data,
                      colors = c("red", "green"),
                      xlab = "Balance Before/After Matching")

# base R graphics are allowed

abline(v = colMeans(balance_data), lty = 3, col = "grey")

