
library(dyn)
set.seed(1)
x <- matrix(rnorm(40), ncol = 2, dimnames = list(NULL,c("gdp", "ipi")))
x[1:10,2]<-NA

x.ts <- ts(x, start = 1950, freq = 4)

x.ts.lm <- dyn$lm(gdp ~ lag(gdp, -2) + lag(ipi, -1), data = x.ts, 
	na.action = na.exclude)

resid(x.ts.lm)
fitted(x.ts.lm)

