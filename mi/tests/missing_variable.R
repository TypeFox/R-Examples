stopifnot(require(mi))

x <- rnorm(10)
x[1] <- NA
y <- missing_variable(x, type = "continuous")

y <- missing_variable(x, type = "irrelevant")
x <- rep(1, 10)
y <- missing_variable(x, type = "fixed")
x <- rep(1:5, each = 2)
y <- missing_variable(x, type = "group")

x[1] <- NA
y <- missing_variable(x, type = "unordered-categorical")
y <- missing_variable(x, type = "ordered-categorical")
y <- missing_variable(x, type = "interval")

x <- rbinom(10, size = 1, prob = 0.5)
x[1] <- NA
y <- missing_variable(x, type = "binary")
y <- missing_variable(x, type = "grouped-binary", strata = rep(c("A", "B"), each = 5))

x <- runif(10)
x[1] <- NA
y <- missing_variable(x, type = "bounded-continuous", lower = 0, upper = 1)
y <- missing_variable(x, type = "positive-continuous")
y <- missing_variable(x, type = "proportion")
x[which.min(x)] <- 0
y <- missing_variable(x, type = "nonnegative-continuous")
y <- missing_variable(x, type = "SC_proportion")
x[which.max(x)] <- 1
y <- missing_variable(x, type = "SC_proportion")

x <- rpois(10, lambda = 5)
x[1] <- NA
y <- missing_variable(x, type = "count")

