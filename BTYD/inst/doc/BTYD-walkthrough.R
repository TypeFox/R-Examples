## ----fig.path="figure/", label="pnbdCalibrationFit", results="hide", echo=FALSE, include=FALSE----
library(knitr)
opts_chunk$set(comment="#");
library(BTYD);
data(cdnowSummary);
est.params <- cdnowSummary$est.params;
cal.cbs <- cdnowSummary$cbs;
cal.fit <- pnbd.PlotFrequencyInCalibration(est.params, cal.cbs, 7)

## ----message=FALSE, tidy=FALSE-------------------------------------------
cdnowElog <- system.file("data/cdnowElog.csv", package = "BTYD")
elog <- dc.ReadLines(cdnowElog, cust.idx = 2,
                     date.idx = 3, sales.idx = 5)
elog[1:3,]

## ----message=FALSE-------------------------------------------------------
elog$date <- as.Date(elog$date, "%Y%m%d");
elog[1:3,]

## ----results="hide", message=FALSE---------------------------------------
elog <- dc.MergeTransactionsOnSameDate(elog);

## ----message=FALSE-------------------------------------------------------
end.of.cal.period <- as.Date("1997-09-30")
elog.cal <- elog[which(elog$date <= end.of.cal.period), ]

## ----results="hide", message=FALSE---------------------------------------
split.data <- dc.SplitUpElogForRepeatTrans(elog.cal);
clean.elog <- split.data$repeat.trans.elog;

## ----message=FALSE-------------------------------------------------------
freq.cbt <- dc.CreateFreqCBT(clean.elog);
freq.cbt[1:3,1:5]

## ----results="hide", message=FALSE---------------------------------------
tot.cbt <- dc.CreateFreqCBT(elog)
cal.cbt <- dc.MergeCustomers(tot.cbt, freq.cbt)

## ----tidy=FALSE, results="hide", message=FALSE---------------------------
birth.periods <- split.data$cust.data$birth.per
last.dates <- split.data$cust.data$last.date
cal.cbs.dates <- data.frame(birth.periods, last.dates, 
                            end.of.cal.period)
cal.cbs <- dc.BuildCBSFromCBTAndDates(cal.cbt, cal.cbs.dates, 
                                      per="week")

## ------------------------------------------------------------------------
params <- pnbd.EstimateParameters(cal.cbs);
params
LL <- pnbd.cbs.LL(params, cal.cbs);
LL

## ------------------------------------------------------------------------
p.matrix <- c(params, LL);
for (i in 1:2){
params <- pnbd.EstimateParameters(cal.cbs, params);
LL <- pnbd.cbs.LL(params, cal.cbs);
p.matrix.row <- c(params, LL);
p.matrix <- rbind(p.matrix, p.matrix.row);
}
colnames(p.matrix) <- c("r", "alpha", "s", "beta", "LL");
rownames(p.matrix) <- 1:3;
p.matrix;

## ----fig.path="figure/", label="pnbdTransactionHeterogeneity", results="hide", include=FALSE----
pnbd.PlotTransactionRateHeterogeneity(params);

## ----fig.path="figure/", label="pnbdDropoutHeterogeneity", results="hide", include=FALSE----
pnbd.PlotDropoutRateHeterogeneity(params);

## ------------------------------------------------------------------------
pnbd.Expectation(params, t=52);

## ----tidy=FALSE----------------------------------------------------------
cal.cbs["1516",]
x <- cal.cbs["1516", "x"]
t.x <- cal.cbs["1516", "t.x"]
T.cal <- cal.cbs["1516", "T.cal"]
pnbd.ConditionalExpectedTransactions(params, T.star = 52, 
                                     x, t.x, T.cal)
pnbd.PAlive(params, x, t.x, T.cal)

## ----tidy=FALSE----------------------------------------------------------
for (i in seq(10, 25, 5)){
  cond.expectation <- pnbd.ConditionalExpectedTransactions(
                        params, T.star = 52, x = i,
                        t.x = 20, T.cal = 39)
  cat ("x:",i,"\t Expectation:",cond.expectation, fill = TRUE)
}

## ----results="hide", eval=FALSE------------------------------------------
#  pnbd.PlotFrequencyInCalibration(params, cal.cbs, 7)

## ----message=FALSE-------------------------------------------------------
elog <- dc.SplitUpElogForRepeatTrans(elog)$repeat.trans.elog;
x.star <- rep(0, nrow(cal.cbs));
cal.cbs <- cbind(cal.cbs, x.star);
elog.custs <- elog$cust;
for (i in 1:nrow(cal.cbs)){
  current.cust <- rownames(cal.cbs)[i]
  tot.cust.trans <- length(which(elog.custs == current.cust))
  cal.trans <- cal.cbs[i, "x"]
  cal.cbs[i, "x.star"] <-  tot.cust.trans - cal.trans
}
cal.cbs[1:3,]

## ----fig.path="figure/", label="pnbdCondExpComp", tidy=FALSE, echo=TRUE, size="small", fig=FALSE----
T.star <- 39 # length of the holdout period
censor <- 7  # This censor serves the same purpose described above
x.star <- cal.cbs[,"x.star"]
comp <- pnbd.PlotFreqVsConditionalExpectedFrequency(params, T.star,
                                              cal.cbs, x.star, censor)
rownames(comp) <- c("act", "exp", "bin")
comp

## ------------------------------------------------------------------------
tot.cbt <- dc.CreateFreqCBT(elog)
d.track.data <- rep(0, 7 * 78)
origin <- as.Date("1997-01-01")
for (i in colnames(tot.cbt)){
  date.index <- difftime(as.Date(i), origin) + 1;
  d.track.data[date.index] <- sum(tot.cbt[,i]);
}
w.track.data <-  rep(0, 78)
for (j in 1:78){
  w.track.data[j] <- sum(d.track.data[(j*7-6):(j*7)])
}

## ----fig.path="figure/", label="pnbdTrackingInc", tidy=FALSE, echo=TRUE, fig=FALSE----
T.cal <- cal.cbs[,"T.cal"]
T.tot <- 78
n.periods.final <- 78
inc.tracking <- pnbd.PlotTrackingInc(params, T.cal, 
                                     T.tot, w.track.data,
                                     n.periods.final)
inc.tracking[,20:25]

## ----fig.path="figure/", label="pnbdTrackingCum", tidy=FALSE, echo=TRUE, fig=FALSE----
cum.tracking.data <- cumsum(w.track.data)
cum.tracking <- pnbd.PlotTrackingCum(params, T.cal, 
                                     T.tot, cum.tracking.data,
                                     n.periods.final)
cum.tracking[,20:25]

## ----fig.path="figure/", label="bgnbdCalibrationFit", results="hide", echo=FALSE, include=FALSE----
data(cdnowSummary);
est.params <- c(0.243, 4.414, 0.793, 2.426);
cal.cbs <- cdnowSummary$cbs;
cal.fit <- bgnbd.PlotFrequencyInCalibration(est.params, cal.cbs, 7)

## ----message=FALSE, tidy=FALSE-------------------------------------------
cdnowElog <- system.file("data/cdnowElog.csv", package = "BTYD")
elog <- dc.ReadLines(cdnowElog, cust.idx = 2,
                     date.idx = 3, sales.idx = 5)
elog[1:3,]

## ----message=FALSE-------------------------------------------------------
elog$date <- as.Date(elog$date, "%Y%m%d");
elog[1:3,]

## ----results="hide", message=FALSE---------------------------------------
elog <- dc.MergeTransactionsOnSameDate(elog);

## ----message=FALSE-------------------------------------------------------
end.of.cal.period <- as.Date("1997-09-30")
elog.cal <- elog[which(elog$date <= end.of.cal.period), ]

## ----results="hide", message=FALSE---------------------------------------
split.data <- dc.SplitUpElogForRepeatTrans(elog.cal);
clean.elog <- split.data$repeat.trans.elog;

## ----message=FALSE-------------------------------------------------------
freq.cbt <- dc.CreateFreqCBT(clean.elog);
freq.cbt[1:3,1:5]

## ----results="hide", message=FALSE---------------------------------------
tot.cbt <- dc.CreateFreqCBT(elog)
cal.cbt <- dc.MergeCustomers(tot.cbt, freq.cbt)

## ----tidy=FALSE, results="hide", message=FALSE---------------------------
birth.periods <- split.data$cust.data$birth.per
last.dates <- split.data$cust.data$last.date
cal.cbs.dates <- data.frame(birth.periods, last.dates, 
                            end.of.cal.period)
cal.cbs <- dc.BuildCBSFromCBTAndDates(cal.cbt, cal.cbs.dates, 
                                      per="week")

## ------------------------------------------------------------------------
params <- bgnbd.EstimateParameters(cal.cbs);
params
LL <- bgnbd.cbs.LL(params, cal.cbs);
LL

## ------------------------------------------------------------------------
p.matrix <- c(params, LL);
for (i in 1:2){
params <- bgnbd.EstimateParameters(cal.cbs, params);
LL <- bgnbd.cbs.LL(params, cal.cbs);
p.matrix.row <- c(params, LL);
p.matrix <- rbind(p.matrix, p.matrix.row);
}
colnames(p.matrix) <- c("r", "alpha", "a", "b", "LL");
rownames(p.matrix) <- 1:3;
p.matrix;

## ----fig.path="figure/", label="bgnbdTransactionHeterogeneity", results="hide", include=FALSE----
bgnbd.PlotTransactionRateHeterogeneity(params);

## ----fig.path="figure/", label="bgnbdDropoutHeterogeneity", results="hide", include=FALSE----
bgnbd.PlotDropoutRateHeterogeneity(params);

## ------------------------------------------------------------------------
bgnbd.Expectation(params, t=52);

## ----tidy=FALSE----------------------------------------------------------
cal.cbs["1516",]
x <- cal.cbs["1516", "x"]
t.x <- cal.cbs["1516", "t.x"]
T.cal <- cal.cbs["1516", "T.cal"]
bgnbd.ConditionalExpectedTransactions(params, T.star = 52, 
                                     x, t.x, T.cal)
bgnbd.PAlive(params, x, t.x, T.cal)

## ----tidy=FALSE----------------------------------------------------------
for (i in seq(10, 25, 5)){
  cond.expectation <- bgnbd.ConditionalExpectedTransactions(
                        params, T.star = 52, x = i,
                        t.x = 20, T.cal = 39)
  cat ("x:",i,"\t Expectation:",cond.expectation, fill = TRUE)
}

## ----results="hide", eval=FALSE------------------------------------------
#  bgnbd.PlotFrequencyInCalibration(params, cal.cbs, 7)

## ----message=FALSE-------------------------------------------------------
elog <- dc.SplitUpElogForRepeatTrans(elog)$repeat.trans.elog;
x.star <- rep(0, nrow(cal.cbs));
cal.cbs <- cbind(cal.cbs, x.star);
elog.custs <- elog$cust;
for (i in 1:nrow(cal.cbs)){
  current.cust <- rownames(cal.cbs)[i]
  tot.cust.trans <- length(which(elog.custs == current.cust))
  cal.trans <- cal.cbs[i, "x"]
  cal.cbs[i, "x.star"] <-  tot.cust.trans - cal.trans
}
cal.cbs[1:3,]

## ----fig.path="figure/", label="bgnbdCondExpComp", tidy=FALSE, echo=TRUE, size="small", fig=FALSE----
T.star <- 39 # length of the holdout period
censor <- 7  # This censor serves the same purpose described above
x.star <- cal.cbs[,"x.star"]
comp <- bgnbd.PlotFreqVsConditionalExpectedFrequency(params, T.star,
                                              cal.cbs, x.star, censor)
rownames(comp) <- c("act", "exp", "bin")
comp

## ------------------------------------------------------------------------
tot.cbt <- dc.CreateFreqCBT(elog)
d.track.data <- rep(0, 7 * 78)
origin <- as.Date("1997-01-01")
for (i in colnames(tot.cbt)){
  date.index <- difftime(as.Date(i), origin) + 1;
  d.track.data[date.index] <- sum(tot.cbt[,i]);
}
w.track.data <-  rep(0, 78)
for (j in 1:78){
  w.track.data[j] <- sum(d.track.data[(j*7-6):(j*7)])
}

## ----results="hide", eval=FALSE------------------------------------------
#  T.cal <- cal.cbs[,"T.cal"]
#  T.tot <- 78
#  n.periods.final <- 78
#  inc.tracking <- bgnbd.PlotTrackingInc(params, T.cal,
#                                       T.tot, w.track.data,
#                                       n.periods.final)
#  inc.tracking[,20:25]

## ----results="hide", eval=FALSE------------------------------------------
#  cum.tracking.data <- cumsum(w.track.data)
#  cum.tracking <- bgnbd.PlotTrackingCum(params, T.cal,
#                                       T.tot, cum.tracking.data,
#                                       n.periods.final)
#  cum.tracking[,20:25]

## ----fig.path="figure/", label="bgbbCalibrationFit", results="hide", echo=FALSE, include=FALSE----
data(donationsSummary);
rf.matrix <- donationsSummary$rf.matrix;
params <- bgbb.EstimateParameters(rf.matrix);
cal.fit <- bgbb.PlotFrequencyInCalibration(params, rf.matrix, 6)

## ----message=FALSE, tidy=FALSE-------------------------------------------
simElog <- system.file("data/discreteSimElog.csv", 
                       package = "BTYD")
elog <- dc.ReadLines(simElog, cust.idx = 1, date.idx = 2)
elog[1:3,]
elog$date <- as.Date(elog$date, "%Y-%m-%d")

max(elog$date);
min(elog$date);
# let's make the calibration period end somewhere in-between
T.cal <- as.Date("1977-01-01")

simData <- dc.ElogToCbsCbt(elog, per="year", T.cal)
cal.cbs <- simData$cal$cbs

freq<- cal.cbs[,"x"]
rec <- cal.cbs[,"t.x"]
trans.opp <- 7 # transaction opportunities
cal.rf.matrix <- dc.MakeRFmatrixCal(freq, rec, trans.opp)
cal.rf.matrix[1:5,]

## ------------------------------------------------------------------------
data(donationsSummary);
rf.matrix <- donationsSummary$rf.matrix
params <- bgbb.EstimateParameters(rf.matrix);
LL <- bgbb.rf.matrix.LL(params, rf.matrix);
p.matrix <- c(params, LL);
for (i in 1:2){
params <- bgbb.EstimateParameters(rf.matrix, params);
LL <- bgbb.rf.matrix.LL(params, rf.matrix);
p.matrix.row <- c(params, LL);
p.matrix <- rbind(p.matrix, p.matrix.row);
}
colnames(p.matrix) <- c("alpha", "beta", "gamma", "delta", "LL");
rownames(p.matrix) <- 1:3;
p.matrix;

## ----fig.path="figure/", label="bgbbTransactionHeterogeneity", results="hide", include=FALSE----
bgbb.PlotTransactionRateHeterogeneity(params);

## ----fig.path="figure/", label="bgbbDropoutHeterogeneity", results="hide", include=FALSE----
bgbb.PlotDropoutRateHeterogeneity(params);

## ------------------------------------------------------------------------
bgbb.Expectation(params, n=10);

## ----tidy=FALSE----------------------------------------------------------
# customer A
n.cal = 6
n.star = 10
x = 0
t.x = 0
bgbb.ConditionalExpectedTransactions(params, n.cal, 
                                     n.star, x, t.x)
# customer B
x = 4
t.x = 5
bgbb.ConditionalExpectedTransactions(params, n.cal, 
                                     n.star, x, t.x)

## ----results="hide", eval=FALSE------------------------------------------
#  bgbb.PlotFrequencyInCalibration(params, rf.matrix)

## ------------------------------------------------------------------------
holdout.cbs <- simData$holdout$cbs
x.star <- holdout.cbs[,"x.star"]

## ----fig.path="figure/", label="bgbbCondExpComp", tidy=FALSE, echo=TRUE, size="small", fig=FALSE----
n.star <- 5 # length of the holdout period
x.star <- donationsSummary$x.star
comp <- bgbb.PlotFreqVsConditionalExpectedFrequency(params, n.star, 
                                                    rf.matrix, x.star)
rownames(comp) <- c("act", "exp", "bin")
comp

## ----fig.path="figure/", label="bgbbCondExpCompRec", tidy=FALSE, echo=TRUE, size="small", fig=FALSE----
comp <- bgbb.PlotRecVsConditionalExpectedFrequency(params, n.star, 
                                                   rf.matrix, x.star)
rownames(comp) <- c("act", "exp", "bin")
comp

## ----fig.path="figure/", label="bgbbTrackingInc", tidy=FALSE, echo=TRUE, fig=FALSE----
inc.track.data <- donationsSummary$annual.trans
n.cal <- 6
xtickmarks <- 1996:2006
inc.tracking <- bgbb.PlotTrackingInc(params, rf.matrix, 
                                     inc.track.data, 
                                     xticklab = xtickmarks)
rownames(inc.tracking) <- c("act", "exp")
inc.tracking

## ----fig.path="figure/", label="bgbbTrackingCum", tidy=FALSE, echo=TRUE, size="small", fig=FALSE----
cum.track.data <- cumsum(inc.track.data)
cum.tracking <- bgbb.PlotTrackingCum(params, rf.matrix, cum.track.data, 
                                     xticklab = xtickmarks)
rownames(cum.tracking) <- c("act", "exp")
cum.tracking

