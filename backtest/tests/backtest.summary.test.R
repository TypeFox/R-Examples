################################################################################
##
## $Id: backtest.summary.test.R 1300 2008-08-27 21:01:11Z zhao $
##
## Tests to make sure the summary method doesn't crash
##
################################################################################

library(backtest)

load("backtest.summary.test.RData")

## save(x, file = "backtest.summary.test.RData", compress = TRUE)

bt.1 <- backtest(x, in.var = "in.var.1", ret.var = "ret.var.1", by.period = FALSE)
bt.2 <- backtest(x, in.var = "in.var.1", ret.var = "ret.var.1", by.var = "country", by.period = FALSE)
bt.3 <- backtest(x, id.var = "id", in.var = "in.var.1", ret.var = "ret.var.1", date.var = "date", by.period = FALSE)
bt.4 <- backtest(x, in.var = c("in.var.1", "in.var.2"), ret.var = "ret.var.1", by.var = "country", by.period = FALSE)
bt.5 <- backtest(x, in.var = c("in.var.1", "in.var.2"), ret.var = c("ret.var.1", "ret.var.2"), by.period = FALSE)
bt.6 <- backtest(x, in.var = c("in.var.1"), ret.var = c("ret.var.1", "ret.var.2"), by.period = FALSE)

summary(bt.1)
summary(bt.2)
summary(bt.3)
summary(bt.4)
summary(bt.5)
summary(bt.6)
