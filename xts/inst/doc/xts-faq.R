### R code from vignette source 'xts-faq.Rnw'

###################################################
### code chunk number 1: preliminaries
###################################################
library("xts")
Sys.setenv(TZ = "GMT")


###################################################
### code chunk number 2: xts-faq.Rnw:86-89 (eval = FALSE)
###################################################
## filenames <- c("a.csv", "b.csv", "c.csv")
## l <- lapply(filenames, read.csv)
## do.call("rbind", l)


###################################################
### code chunk number 3: xts-faq.Rnw:105-106 (eval = FALSE)
###################################################
## lm(myxts[, "Res"] ~ myxts[, "ThisVar"] + myxts[, "ThatVar"])


###################################################
### code chunk number 4: xts-faq.Rnw:109-110 (eval = FALSE)
###################################################
## with(myxts, lm(Res ~ ThisVar + ThatVar))


###################################################
### code chunk number 5: xts-faq.Rnw:117-121
###################################################
x <- .xts(c(1, 2, 3, 0, 0, 0), 1:6)
x[x==0] <- NA
na.locf(x)
x


###################################################
### code chunk number 6: xts-faq.Rnw:128-130
###################################################
data(sample_matrix)
sample.xts = xts(sample_matrix, Sys.time() + seq(0, by = 0.1, length = 180))


###################################################
### code chunk number 7: xts-faq.Rnw:138-140
###################################################
options(digits.secs = 3)
head(sample.xts)


###################################################
### code chunk number 8: xts-faq.Rnw:151-152
###################################################
print(as.numeric(as.POSIXlt("2012-03-20 09:02:50.001")), digits = 20)


###################################################
### code chunk number 9: xts-faq.Rnw:161-162 (eval = FALSE)
###################################################
## myxts.2 <- xts(t(apply(myxts, 1 , myfun)), index(myxts))


###################################################
### code chunk number 10: xts-faq.Rnw:173-174 (eval = FALSE)
###################################################
## do.call(rbind, lapply(split(myxts,"days"), myfun))


###################################################
### code chunk number 11: xts-faq.Rnw:180-182 (eval = FALSE)
###################################################
## rt <- r['T16:00/T17:00','Value']
## rd <- apply.daily(rt, function(x) xts(t(quantile(x,0.9)), end(x)))


###################################################
### code chunk number 12: xts-faq.Rnw:191-199 (eval = FALSE)
###################################################
## # align index into 3-hour blocks
## a <- align.time(s, n=60*60*3)
## # find the number of obs in each block
## count <- period.apply(a, endpoints(a, "hours", 3), length)
## # create an empty \pkg{xts} object with the desired index
## e <- xts(,seq(start(a),end(a),by="3 hours"))
## # merge the counts with the empty object and fill with zeros
## out <- merge(e,count,fill=0)


###################################################
### code chunk number 13: xts-faq.Rnw:209-210 (eval = FALSE)
###################################################
## myxts = as.xts(transform(myxts, ABC = 1))


###################################################
### code chunk number 14: xts-faq.Rnw:214-215 (eval = FALSE)
###################################################
## indexTZ(myxts) = Sys.getenv("TZ")


###################################################
### code chunk number 15: xts-faq.Rnw:227-228 (eval = FALSE)
###################################################
## myts[myts$Symbol == "AAPL" & index(myts) == as.POSIXct("2011-09-21"),]


###################################################
### code chunk number 16: xts-faq.Rnw:231-232 (eval = FALSE)
###################################################
## myts[myts$Symbol == "AAPL"]['2011-09-21']


###################################################
### code chunk number 17: xts-faq.Rnw:239-243
###################################################
data(sample_matrix)
sample.xts <- as.xts(sample_matrix, descr='my new xts object')
x <- sample.xts['2007']  
x[.indexwday(x) %in% 1:5]


###################################################
### code chunk number 18: xts-faq.Rnw:254-255 (eval = FALSE)
###################################################
## qxts = xts(q[,-1], order.by=q[,1])


###################################################
### code chunk number 19: xts-faq.Rnw:268-271 (eval = FALSE)
###################################################
## xTemps <- align.time(xts(temps[,2],as.POSIXct(temps[,1])), n=600)
## xGas <- align.time(xts(gas[,2],as.POSIXct(gas[,1])), n=600)
## merge(xTemps,xGas)


