
library(dyn)

cat("# create some data\n")
set.seed(1)
x <- matrix(rnorm(40), ncol = 2, dimnames = list(NULL, c("gdp", "ipi")))
x.ts <- ts(x, start = 1950, freq = 4)


cat("Run a dynamic regression\n")
x.ts.lm <- dyn$lm(gdp ~ lag(gdp, -2) + lag(ipi,-1), data = x.ts)

resid(x.ts.lm)
fitted(x.ts.lm)

x.pred <- predict(x.ts.lm, x.ts)
x.pred

cat("Plot gdp in black and predicted gdp in red\n")
ts.plot( x.ts[,"gdp"], x.pred, 
	gpars = list(type = "b", col = 1:2, pch = 20))


cat("Run demo() to find names of additional dyn demos\n")

