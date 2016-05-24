# pmcm.R -- version 2010-12-05

x <- rnorm(20L)

# mean and probability of loss
theta     <- 0.1 # ... or mean(x)
prob.loss <- ecdf(x)(theta)
exponent  <- 2

# conditional moment (CM)
(cm <- mean((x[x < theta] - theta)^exponent))

# partial moment (PM)
xx <- x - theta; xx[xx > 0] <- 0
(pm <- mean(xx^exponent))

# relationship between PM and CM
stopifnot(all.equal(cm * prob.loss, pm))