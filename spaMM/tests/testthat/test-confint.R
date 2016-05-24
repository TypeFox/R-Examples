# confint
cat("\ntest confint:")
data(wafers)
wfit <- HLfit(y ~X1+(1|batch),family=Gamma(log),data=wafers,HLmethod="ML")
ci <- confint(wfit,"X1")

expect_equal(ci$interval[[1]],0.01828659,tolerance=1e-4)
expect_equal(ci$interval[[2]],0.17271333,tolerance=1e-4)