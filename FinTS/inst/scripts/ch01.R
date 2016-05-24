### ch. 1.  Financial Time Series and Their Characteristics
###
### 
### Ruey S. Tsay (2005)
### Analysis of Financial Time Series, 2nd ed.
### (Wiley)
###
### 
library(FinTS)

##
## sec. 1.1.  Assett Returns
##

# p. 4
# Table 1.1.  Illustration of the effects of compounding
freqs <- c(1, 2, 4, 12, 52, 365, Inf)
cI <- compoundInterest(0.1, frequency=freqs)
(Table1.1 <- data.frame(Type=c(
    "Annual", "Semiannual", "Quarterly",
    "Monthly", "Weekly", "Daily", "Continuously"),
  Number.of.payments=freqs, 
  Interest.rate.per.period=0.1/freqs,
  Net.Value=cI) )

# p.  6  
# Example 1.1.
logRtn <- .0446
100*expm1(logRtn)

mo.logRtns <- c(.0446, -.0734, .1077)
(q.logRtns <- sum(mo.logRtns))

##
## sec. 1.2.  Distributional Properties of Returns
##
# p. 7
# sec. 1.2.1.  Review of Statistical Distributions and their Moments

# p. 10
# Example 1.2.
data(d.ibmvwewsp6203)
s.ibm <- sd(d.ibmvwewsp6203[, "IBM"])
(skew.ibm <- skewness(d.ibmvwewsp6203[, "IBM"]))
(n.ibm <- dim(d.ibmvwewsp6203)[1])
(skew.sd <- sqrt(6/n.ibm))
(t.sk <- (skew.ibm/skew.sd))
pnorm(-t.sk)

# A slightly more accurate version of this is
# agostino.test 
library(moments)
(skew.test <- agostino.test(as.numeric(d.ibmvwewsp6203[, "IBM"])))

# p. 11
# Table 1.2.  Descriptive Statistics for Daily and Monthly
#             Simple and Log Returns for Selected Indexes and Stocks
data(d.ibmvwewsp6203)
data(d.intc7303)
data(d.3m6203)
data(d.msft8603)
data(d.c8603)
(Daily.Simple.Returns.pct <- rbind(
    FinTS.stats(100*d.ibmvwewsp6203[, "SP"]),     
    FinTS.stats(100*d.ibmvwewsp6203[, "VW"]),     
    FinTS.stats(100*d.ibmvwewsp6203[, "EW"]), 
    FinTS.stats(100*d.ibmvwewsp6203[, "IBM"]),     
    FinTS.stats(100*d.intc7303[,"Intel"]),     
    FinTS.stats(100*d.3m6203[, "MMM"]), 
    FinTS.stats(100*d.msft8603[, 'MSFT']), 
    FinTS.stats(100*d.c8603[, "C"]) 
) )

(Daily.log.Returns.pct <- rbind(
    FinTS.stats(100*log(1+d.ibmvwewsp6203[, "SP"])), 
    FinTS.stats(100*log(1+d.ibmvwewsp6203[, "VW"])),     
    FinTS.stats(100*log(1+d.ibmvwewsp6203[, "EW"])), 
    FinTS.stats(100*log(1+d.ibmvwewsp6203[, "IBM"])),     
    FinTS.stats(100*log(1+d.intc7303[,"Intel"])),     
    FinTS.stats(100*log(1+d.3m6203[, "MMM"])), 
    FinTS.stats(100*log(1+d.msft8603[, 'MSFT'])), 
    FinTS.stats(100*log(1+d.c8603[, "C"])) 
) )

data(m.ibmvwewsp2603)
data(m.intc7303)
data(m.3m4603)
data(m.msft8603)
data(m.c8603)
(Monthly.Simple.Returns.pct <- rbind(
    SP=FinTS.stats(100*m.ibmvwewsp2603[, "SP"]),     
    VW=FinTS.stats(100*m.ibmvwewsp2603[, "VW"]),     
    EW=FinTS.stats(100*m.ibmvwewsp2603[, "EW"]), 
    IBM=FinTS.stats(100*m.ibmvwewsp2603[, "IBM"]),     
    Intel=FinTS.stats(100*m.intc7303[,"Intel"]),     
    MMM=FinTS.stats(100*m.3m4603[, "MMM"]), 
    Microsoft=FinTS.stats(100*m.msft8603[, 'MSFT']), 
    CitiGroup=FinTS.stats(100*m.c8603[, "C"]) 
) )

(Monthly.log.Returns.pct <- rbind(
    SP=FinTS.stats(100*log(1+m.ibmvwewsp2603[, "SP"])), 
    VW=FinTS.stats(100*log(1+m.ibmvwewsp2603[, "VW"])),     
    EW=FinTS.stats(100*log(1+m.ibmvwewsp2603[, "EW"])), 
    IBM=FinTS.stats(100*log(1+m.ibmvwewsp2603[, "IBM"])),     
    Intel=FinTS.stats(100*log(1+m.intc7303[,"Intel"])),     
    MMM=FinTS.stats(100*log(1+m.3m4603[, "MMM"])), 
    Microsoft=FinTS.stats(100*log(1+m.msft8603[, 'MSFT'])), 
    CitiGroup=FinTS.stats(100*log(1+m.c8603[, "C"])) 
) )

dimnames(Daily.Simple.Returns.pct)[[1]] <-
  c("SP", "VW", "EW", "IBM", "Intel", "3M", "MSFT", "Citi")
dimnames(Daily.log.Returns.pct)[[1]] <-
  c("SP", "VW", "EW", "IBM", "Intel", "3M", "MSFT", "Citi")
dimnames(Monthly.Simple.Returns.pct)[[1]] <-
  c("SP", "VW", "EW", "IBM", "Intel", "3M", "MSFT", "Citi")
dimnames(Monthly.log.Returns.pct)[[1]] <-
  c("SP", "VW", "EW", "IBM", "Intel", "3M", "MSFT", "Citi")

#write.table(rbind(Daily.Simple.Returns.pct, Daily.log.Returns.pct,
#  Monthly.Simple.Returns.pct, Monthly.log.Returns.pct),"Table1-2.csv", sep=",")

# p. 12

# Comparable to S-Plus Demonstration
# module(finmetrics)
# x = matrix(scan(file = 'd-ibmvwewsp6203.txt'), 5) % load the data
# ibm = x[2,]*100 % compute percentage returns

data(d.ibmvwewsp6203)
quantile(100*as.numeric(d.ibmvwewsp6203[, "IBM"]))
FinTS.stats(100*d.ibmvwewsp6203[, "IBM"])

# p. 15
# Sec. 1.2.2.  Distributions of Returns
# p. 16
# Figure 1.1.  Comparisons of Finite Mixture, Stable
#              and Standard Normal Distributions
library(distrEx)

N01 <- Norm()
N04 <- Norm(0, 4)

contamNorm <- ConvexContamination(N01, N04, size = 0.05)

plot(dnorm, xlim=c(-4, 4))
text(.75, .39, "Normal") 
x <- seq(-4, 4, len=201)
lines(x, d(contamNorm)(x), lty="dotted", col="red", lwd=2)
text(0, 0.351, "Mixture", col="red")
lines(x, dcauchy(x), lty="dashed", col="green", lwd=2)
text(0, 0.25, "Cauchy", col="green")

# p. 17
# Sec. 1.2.5.  Empirical properties of returns

# p. 18 
# Figure 1.2.  Time plots of monthly returns of IBM stock

op <- par(mfrow=c(2,1))
plot(m.ibmvwewsp2603[, "IBM"], xlab="year", ylab="s-rtn")
plot(log(1+m.ibmvwewsp2603[, "IBM"]), xlab="year", ylab="log-rtn")
par(op)


# Figure 1.3.  Time plots of monthly returens of the value-weighted index

op <- par(mfrow=c(2,1))
plot(m.ibmvwewsp2603[, "VW"], xlab="year", ylab="s-rtn")
plot(log(1+m.ibmvwewsp2603[, "VW"]), xlab="year", ylab="log-rtn")
par(op)

# p. 19
# Figure 1.4.  Comparison of empirical and normal densities
# for the monthly simple and log returns of IBM stock

if(require(logspline)){

  m.ibm <- m.ibmvwewsp2603[, "IBM"]
  dens.ibm.rtns <- logspline(100*as.numeric(m.ibm))
  m.log.ibm <- 100*log(1+m.ibm)
  dens.ibm.logrtns <- logspline(as.numeric(m.log.ibm))

  op <- par(mfrow=c(1,2))
  plot(dens.ibm.rtns, xlim=c(-40, 40), xlab="simple returns", ylab="density")
  x.ibm <- seq(-40, 40, len=201)
  lines(x.ibm, dnorm(x.ibm, mean=100*mean(m.ibm), s=100*sd(m.ibm)),
        lty="dotted", col="red", lwd=2)

  plot(dens.ibm.logrtns, xlim=c(-40, 40), xlab="log returns", ylab="density")
  lines(x.ibm, dnorm(x.ibm, mean=mean(m.log.ibm), s=sd(m.log.ibm)),
        lty="dotted", col="red", lwd=2)
  par(op)
}
qqnorm(100*m.ibm, datax=TRUE)
qqnorm(m.log.ibm, datax=TRUE)

# normal plots may provide a more sensitive evaluation
# of skewness and kurtosis than density plots

# The kurtosis here does not look extreme.
# Similar plots of the daily returns show
# much more extreme kurtosis.

# p. 20
##
## 1.3.  Processes Considered
##
data(m.gs10)
data(m.gs1)
# Figure 1.5.  Time plot of US monthly interest rates 
op <- par(mfrow=c(2,1))
plot(m.gs10, xlab="year", ylab="rate", type="l")
title("Monthly 10-yr Treasury constant maturity rates")
plot(m.gs1, xlab="year", ylab="rate", type="l")
title("Monthly 1-yr Treasury constant maturity rates")
par(op)

# p. 21
# Figure 1.6.  Time plot of daily exchange rate
#              between US dollara and Japanese Yen 
data(d.fxjp00)
op <- par(mfrow=c(2,1))
plot(d.fxjp00, xlab="year", ylab="Yens", type="l")
plot(diff(d.fxjp00), xlab="year", ylab="Change", type="l")
par(op)

# p. 22
# Table 1.3.  Descriptive Statistics of Selected US Financial Time Series
data(m.fama.bond5203)
(m.bondRtns <- rbind(
  "1-12 months"=FinTS.stats(100*m.fama.bond5203[, "m1.12"]), 
  "24-36 months"=FinTS.stats(100*m.fama.bond5203[, "m24.36"]), 
  "48-60 months"=FinTS.stats(100*m.fama.bond5203[, "m48.60"]), 
  "61-120 months"=FinTS.stats(100*m.fama.bond5203[, "m61.120"]) ))

data(m.gs1)
data(m.gs3)
data(m.gs5)
data(m.gs10)
(m.treasuryRtns <- rbind(
  "1 year"=FinTS.stats(m.gs1),
  "3 years"=FinTS.stats(m.gs3),
  "5 years"=FinTS.stats(m.gs5),
  "10 years"=FinTS.stats(m.gs10) ))

data(w.tb3ms)
data(w.tb6ms)
(w.treasuryRtns <- rbind(
  "3 months"=FinTS.stats(w.tb3ms),
  "6 months"=FinTS.stats(w.tb6ms) ) )

#write.table(rbind(m.bondRtns, m.treasuryRtns, w.treasuryRtns),
#            "Table1-3.csv", sep=",")
