### R code from vignette source 'backtest.Rnw'

###################################################
### code chunk number 1: backtest.Rnw:24-36
###################################################
options(width = 50, digits = 2, scipen = 5)
set.seed(1)
cat.df.without.rownames <- function (d, file = ""){
  stopifnot(is.data.frame(d))
  row.names(d) <- 1:nrow(d)
  x <- NULL
  conn <- textConnection("x", "w", local = TRUE)
  capture.output(print(d), file = conn)
  close(conn)
  cat(substring(x, first = max(nchar(row.names(d))) + 2), sep = "\n", 
      file = file)
}


###################################################
### code chunk number 2: backtest.Rnw:96-97
###################################################
library(backtest)


###################################################
### code chunk number 3: backtest.Rnw:100-102
###################################################
data(starmine)
names(starmine)


###################################################
### code chunk number 4: backtest.Rnw:112-113
###################################################
cat.df.without.rownames(as.data.frame.table(table(date = as.character(starmine$date)), responseName = "count"))


###################################################
### code chunk number 5: backtest.Rnw:120-121
###################################################
current.options <- options(digits = 1, width = 80, scipen = 99)


###################################################
### code chunk number 6: backtest.Rnw:124-125
###################################################
cat.df.without.rownames(starmine[row.names(starmine) %in% c(2254, 9852, 10604, 18953, 62339, 77387, 85739), c("date", "name", "ret.0.1.m", "ret.0.6.m","smi")])


###################################################
### code chunk number 7: backtest.Rnw:128-129
###################################################
options(current.options)


###################################################
### code chunk number 8: backtest.Rnw:143-144
###################################################
bt <- backtest(starmine, in.var = "smi", ret.var = "ret.0.1.m", by.period = FALSE)


###################################################
### code chunk number 9: backtest.Rnw:157-161
###################################################
table(cut(starmine$smi, breaks = quantile(starmine$smi, probs=seq(0,1,0.20), 
                          na.rm=TRUE, 
                          names=TRUE), 
          include.lowest = TRUE))


###################################################
### code chunk number 10: backtest.Rnw:173-174
###################################################
summary(bt)


###################################################
### code chunk number 11: natural portfolio
###################################################
bt <- backtest(starmine, id.var = "id",
               date.var = "date", in.var = "smi",
               ret.var = "ret.0.1.m", natural = TRUE, by.period = FALSE)


###################################################
### code chunk number 12: backtest.Rnw:232-234
###################################################
current.options <- options()
options(digits = 1, width = 80)


###################################################
### code chunk number 13: summary natural backtest
###################################################
summary(bt)


###################################################
### code chunk number 14: backtest.Rnw:245-246
###################################################
options(current.options)


###################################################
### code chunk number 15: backtest.Rnw:302-303
###################################################
op <- options(digits = 1)


###################################################
### code chunk number 16: smi.2 backtest
###################################################
bt <- backtest(starmine,
                     id.var = "id",
                     date.var = "date",
                     in.var = c("smi", "cap.usd"),
                     ret.var = "ret.0.1.m",
                     natural = TRUE,
                     by.period = FALSE)


###################################################
### code chunk number 17: backtest.Rnw:316-317
###################################################
bt.save <- bt


###################################################
### code chunk number 18: backtest.Rnw:324-325
###################################################
summary(bt)


###################################################
### code chunk number 19: backtest.Rnw:327-328
###################################################
options(op)


###################################################
### code chunk number 20: backtest.Rnw:348-349
###################################################
plot(bt, type = "return")


###################################################
### code chunk number 21: backtest.Rnw:363-364
###################################################
plot(bt.save, type = "cumreturn.split")


###################################################
### code chunk number 22: backtest.Rnw:385-386
###################################################
plot(bt, type = "turnover")


###################################################
### code chunk number 23: backtest.Rnw:404-406
###################################################
current.options <- options()
options(digits = 1, width = 50)


###################################################
### code chunk number 24: sector backtest
###################################################
bt <- backtest(starmine, in.var = "smi", ret.var = "ret.0.1.m", by.var = "sector", by.period = FALSE)


###################################################
### code chunk number 25: backtest.Rnw:413-414
###################################################
options(width = 80)


###################################################
### code chunk number 26: backtest.Rnw:417-418
###################################################
summary(bt)


###################################################
### code chunk number 27: backtest.Rnw:421-422
###################################################
options(current.options)


###################################################
### code chunk number 28: counts
###################################################
counts(bt)


###################################################
### code chunk number 29: backtest.Rnw:452-454
###################################################
current.options <- options()
options(digits = 1, width = 30)


###################################################
### code chunk number 30: market cap
###################################################
bt <- backtest(starmine, in.var = "smi", ret.var = "ret.0.1.m", by.var = "cap.usd", buckets = c(5, 10), by.period = FALSE)


###################################################
### code chunk number 31: backtest.Rnw:461-462
###################################################
options(current.options)


###################################################
### code chunk number 32: backtest.Rnw:465-466
###################################################
summary(bt)


###################################################
### code chunk number 33: backtest.Rnw:486-488
###################################################
bt <- backtest(starmine, in.var = "smi", buckets = 4,
               ret.var = c("ret.0.1.m", "ret.0.6.m"), by.period = FALSE)


###################################################
### code chunk number 34: backtest.Rnw:491-492
###################################################
op <- options(width = 80)


###################################################
### code chunk number 35: backtest.Rnw:495-496
###################################################
summary(bt)


###################################################
### code chunk number 36: backtest.Rnw:511-512
###################################################
options(op)


