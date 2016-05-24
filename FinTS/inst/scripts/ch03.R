### ch. 3.  Conditional Heteroscedastic Models 
###
### 
### Ruey S. Tsay (2005)
### Analysis of Financial Time Series, 2nd ed.
### (Wiley)
###
### 

# p. 97
library(FinTS)

# p. 98 
##
## sec. 3.1.  Characteristics of volatility 
##

# p. 99
##
## sec. 3.2.  Structure of a model 
##

# Figure 3.1
data(m.intc7303)
ml.intc <- log(1+m.intc7303)

op <- par(mfcol=c(2,2))
Acf(as.numeric(ml.intc), main="(a) Log returns")
Acf(as.numeric(ml.intc)^2, main="(b) Squared log returns")
Acf(abs(as.numeric(ml.intc)), main="(c) Absolute log returns") 
pacf(as.numeric(ml.intc)^2, main="(d) Squared log returns")
par(op)
# adjust ylim = c(-.2, .4)
op <- par(mfcol=c(2,2))
Acf(as.numeric(ml.intc), main="(a) Log returns", ylim = c(-.2, .4))
Acf(as.numeric(ml.intc)^2, main="(b) Squared log returns", ylim = c(-.2, .4))
Acf(abs(as.numeric(ml.intc)), main="(c) Absolute log returns", ylim = c(-.2, .4)) 
pacf(as.numeric(ml.intc)^2, main="(d) Squared log returns", ylim = c(-.2, .4))
par(op) 


# p. 101 
##
## sec. 3.3.  Model Building 
##

# sec. 3.3.1.  Testing for ARCH effect 

# p. 102
data(m.intc7303)
str(m.intc7303)

LjB.intc <- AutocorTest(log(1+as.numeric(m.intc7303)), lag=12)

LM.intc <- ArchTest(log(1+as.numeric(m.intc7303)), lag=12)

##
## sec. 3.4.  The ARCH model
##

# p. 103 
# Figure 3.2
data(exch.perc)
str(exch.perc)

op <- par(mfrow=c(2,1))
plot(exch.perc, type="l", xlab="", ylab="fx",
     main="(a) Percentage change in exchange rate")
plot(exch.perc^2, type="l", xlab="", ylab="sq-fx",
     main="(b) Squared series")
par(op)

# p. 104 
# Figure 3.3

op <- par(mfrow=c(2,1))
Acf(exch.perc, ylim=c(-.1, .1), main="(a) Sample ACF")

pacf(exch.perc^2, ylim=c(-.1, .1),
     main="(b) Partial ACF of the squared series")
par(op)

# sec. 3.4.1.  Properties of ARCH models   

# p. 106
# sec. 3.4.2.  Weaknesses of ARCH models

# sec. 3.4.3.  Building an ARCH model 

# p. 109
# sec. 3.4.4.  Some Examples

# Example 3.1
data(m.intc7303)
ml.intc <- log(1+m.intc7303)

# Possibiliities:
#garch {tseries} Fit GARCH Models to Time Series
#GarchFitting {fGarch} Univariate GARCH Time Series Fitting

library(tseries)
#arch3.fit <- garch(ml.intc, order=c(0, 3))
#Error in garch(ml.intc, order = c(1, 3)) : NAs in x
# introduced by conversion from class 'zoo'

(ml.intc. <- mean(ml.intc))
# Moderatly close to "C" in the garch(3, 0) and garch(1, 0) results on p. 103
ml.intc0 <- ml.intc-mean(ml.intc)
arch3.fit <- garch(as.numeric(ml.intc0), order=c(0, 3))
summary(arch3.fit)
# Not a great match, but moderately close
# to the garch(3, 0) results on p. 103.

(arch1.fit <- garch(as.numeric(ml.intc0), order=c(0, 1)))
# Closer to the book than for arch3.fit

stdresi <- residuals(arch1.fit)

AutocorTest(stdresi, 10)

ArchTest(stdresi, 10)

# unconditional standard deviation:  
sqrt(with(arch1.fit, coef[1]/(1-coef[2])))


######################################
##
## NEED:  a function to fit ARCH(1) with Student's t innovations (p. 111)  
##
## ... with asypmtotic stdev
##
######################################

library(fGarch)
#arch3.Fit <- garchFit(~garch(3, 0), data=ml.intc)
#Error in sum(beta) : invalid 'type' (closure) of argument

##??????????????
##
## Questions sent to a package maintainer.
## 2007.12.24 and 2008.02.16 
##
##??????????????

resi <- residuals(arch1.fit)

AutocorTest(resi[-1], 10)
AutocorTest(resi[-1]^2, 10)

##############################
##
## t Innovations?
##
## fitARCH.t1 <- garchFit(~garch(1, 0), data=ml.intc, cond.dist = 'std')
## Error in sum(beta) : invalid 'type' (closure) of argument ??? 
##
##############################
# p. 112
# Figure 3.4

op <- par(mfcol=c(2,2))
Acf(stdresi[-1], main='(a) ARCH(1) standardized residuals')
Acf(stdresi[-1]^2, main='(b) ARCH(1) squared standardized residuals')
Acf(abs(stdresi[-1]), main='(c) ARCH(1) abs(standardized residuals)')
plot(stdresi, type='l', main='(d) standardized residuals')
par(op)

op <- par(mfcol=c(2,2))
Acf(stdresi[-1], main='(a) ARCH(1) standardized residuals', ylim=c(-.2, .4))
Acf(stdresi[-1]^2, main='(b) ARCH(1) squared standardized residuals',
    ylim=c(-.2, .4))
Acf(abs(stdresi[-1]), main='(c) ARCH(1) abs(standardized residuals)'
    , ylim=c(-.2, .4))
plot(stdresi, type='l', main='(d) standardized residuals')
par(op)

# p. 113

# Example 3.2.  
data(exch.perc)

library(tseries)
exch.arch1 <- garch(exch.perc, order=c(0, 3))
summary(exch.arch1)

##
## sec. 3.5.  The GARCH Model
##

# p. 116
#  sec. 3.5.1.  An Illustrative Example
data(sp500)

(spFit03 <- arima(sp500, c(0, 0, 3)))
coef(spFit03)
# not great:  Force ma2 = 0 
(spFit03. <- arima(sp500, c(0, 0, 3),
                 fixed=c(ma1=NA, ma2=0, ma3=NA, intercept=NA)))
names(spFit03.)
sqrt(spFit03.$sigma2)

str(sp500)
(spFit30 <- arima(sp500, c(3, 0, 0)))

library(fGarch)
spFit30.11 <- garchFit(sp500~arma(3,0)+garch(1,1),
                       data=sp500)
spFit30.11
# Difference from the book could be minor differences in
# the accuracy of the nonlinear optimizer used.  

# p. 117

# Figure 3.5
plot(sp500, xlab="year", ylab="rtn")
abline(h=0, lty="dashed")

# Figure 3.6 
op <- par(mfrow=c(2,1))
Acf(sp500, lag.max=30, main="(a)")
pacf(sp500^2, main="(b)")
par(op)

# p. 118
# unconditional var(a[t]) =
#str(spFit30.11)
with(spFit30.11@fit, par["omega"]/(1-par["beta1"]-par["alpha1"]))
# Differs from the book by 10%
# ... within the numeric precision of available
# algorithms for this type of problem?  

# Drop insignificant terms and refit:  

data(sp500)
spFit00.11 <- garchFit(sp500~garch(1,1), data=sp500)
#Warning messages:
#1: In if (class(data) == "timeSeries") { :
#  the condition has length > 1 and only the first element will be used
#2: In if (class(data) == "data.frame") { :
#  the condition has length > 1 and only the first element will be used
                                         
spFit00.11a <- garchFit(sp500~garch(1,1), data=as.numeric(sp500))
# No warnings ... 
spFit00.11
spFit00.11a
all.equal(spFit00.11, spFit00.11a)
#  Answers seem virtually identical 

# unconditional var(a[t]) =
with(spFit00.11@fit, par["omega"]/(1-par["beta1"]-par["alpha1"]))
# again within 10%

str(spFit00.11)

# Figure 3.7 
op <- par(mfrow=c(2,1))
plot(index(sp500), spFit00.11@sigma.t, type="l", xlab="year",
     ylab="sigma.t", main="(a) Estimated volatility process")
std.res <- residuals(spFit00.11)/spFit00.11@sigma.t
plot(index(sp500), std.res, type="l", xlab="year",
     ylab="std-resi", main="(b) Standardized residuals")
par(op)

# p. 119
# Figure 3.8
op <- par(mfrow=c(2,1))
Acf(std.res, ylim=c(-.2, .2), main="(a)", lag.max=24)
Acf(std.res^2, ylim=c(-.2, .2), main="(b)", lag.max=24)
par(op)

AutocorTest(std.res, 12)
AutocorTest(std.res, 24)

AutocorTest(std.res^2, 12)
AutocorTest(std.res^2, 24)

# p. 120
# Since the 'arma' part of the model is nothing,
# the "mean return" forecast is constant 'mu',
# estimated here as 0.0074497 vs. 0.0076 in the book.  
#str(spFit00.11)
names(spFit00.11@fit)
spFit00.11@fit$coef['mu']

pred.spFit00.11 <- predict(spFit00.11, 1000)

plot(pred.spFit00.11[, "standardDeviation"])
pred.spFit00.11[c(1:5, 1000), -2]

# Moderately close to Table 3.1  

# With Student's t innovations, 5 degrees of freedom:  
spFit00.11t5 <- garchFit(sp500~garch(1,1), data=as.numeric(sp500),
                        cond.dist = "std", shape=5, include.shape=FALSE)
spFit00.11t5

#######################################
##
## The answers here seem identical to spFit00.11a
## Problem reported to maintainer 2008.02.23
##
## Code retained in the hopes that 'garchFit' will be fixed
## and give better answers in the future
##
#######################################

str(spFit00.11t5)
sum(spFit00.11t5@fit$coef[c("alpha1", "beta1")])
# 0.9763 ... nearly IGARCH

stdres.ti5 <- (spFit00.11t5@residuals / spFit00.11t5@sigma.t)

AutocorTest(stdres.ti5, lag=10)
AutocorTest(stdres.ti5^2, lag=10)

# estimate the degrees of freedom:  

spFit00.11t <- garchFit(sp500~garch(1,1), data=as.numeric(sp500),
                        cond.dist = "std")
spFit00.11t
spFit00.11t5
# answers roughly within 10% of the book.  

str(spFit00.11t)
sum(spFit00.11t@fit$coef[c("alpha1", "beta1")])
# 0.97633 ... nearly IGARCH

stdres.ti <- (spFit00.11t@residuals / spFit00.11t@sigma.t)

AutocorTest(stdres.ti, lag=10)
AutocorTest(stdres.ti^2, lag=10)


#######################################
##
## The answers here seem identical to spFit00.11a
## Problem reported to maintainer 2008.02.23
##
## Code retained in the hopes that 'garchFit' will be fixed
## and give better answers in the future
##
#######################################

# p. 121
# sec. 3.5.2.  Forecasting evaluation

# sec. 3.5.3.  A two-pass estimation method

mean(sp500)
# 0.006143

sp500a <- sp500-mean(sp500)

(sp500a1.1 <- arima(sp500a^2, c(1,0,1)))

coef(sp500a1.1)

(beta.1 <- (-coef(sp500a1.1)[2]))
(alpha.1 <- sum(coef(sp500a1.1)[1:2]))

# p. 122
##
## 3.6.  The Integrated GARCH model
##

# How to estimate IGARCH?

#############################
##
## fGarch includes 'GarchOxInterface',
## which seems to require the Ox commercial software (???)
##
#############################

# p. 123
##
## 3.7.  The GARCH-M model
##

#############################
##
## How to estimate GARCH-M?
##
## Question sent to R-SIG-FINANCE 2008.02.23
##
#############################

# GARCH(1,1)-M example for SP500 

# p. 124
##
## 3.8.  The exponential GARCH model
##

#############################
##
## How to estimate EGARCH?
##
#############################

# p. 130
##
## 3.9.  Threshold GARCH
##

#############################
##
## How to estimate TGARCH?
##
#############################

# p. 131
##
## 3.10.  CHARMA 
##

#############################
##
## How to estimate CHARMA?
##
#############################

# p. 133
##
## 3.11.  Random Coefficient Autoregressive Models  
##

#############################
##
## How to estimate ?  
##
#############################

# p. 134
##
## 3.12.  Stochastic Volatility
##

#############################
##
## How to estimate ?  
##
#############################

##
## 3.13.  Long-Memory Stochastic Volatility
##

data(d.sp8099)
data(d.ibmvwewsp6203)

# Figure 3.9
op <- par(mfrow=c(2,1))
Acf(abs(log(1+as.numeric(d.sp8099))), 200)
# Similar but different from Figure 3.9(a)
Acf(abs(log(1+as.numeric(d.ibmvwewsp6203[, "IBM"]))), 200)
par(op)

#############################
##
## How to estimate ?  
##
#############################

# p. 136
##
## 3.14.  Applications
##

data(m.ibmvwewsp2603)

op <- par(mfrow=c(2,1))
plot(log(1+m.ibmvwewsp2603[, "IBM"]), type="l",
     xlab="year", ylab="log-rtn", main="(a) IBM")
plot(log(1+m.ibmvwewsp2603[, "SP"]), type="l",
     xlab="year", ylab="log-rtn", main="(a) S & P 500 Index")
par(op)

ibmFit10.11 <- garchFit(~arma(1,0)+garch(1,1),
                     data=as.numeric(m.ibmvwewsp2603[, "IBM"]))
ibmFit10.11
ibmStd.res <- residuals(ibmFit10.11)/ibmFit10.11@sigma.t

AutocorTest(ibmStd.res, 10)
AutocorTest(ibmStd.res, 20)
AutocorTest(ibmStd.res^2, 10)
AutocorTest(ibmStd.res^2, 20)

m2603Dates <- index(m.ibmvwewsp2603)
str(m2603Dates)
# yearmon = year + (i-1)/12
# i = 1:12

# summer
library(seas)

summer2603 <- (mkseas(as.Date(m2603Dates), width='DJF')=='JJA')

################################
##
## How to estimate GARCH with an explanatory variable
##      like summer2603 in the garch model?  
##
################################



# p. 140
##
## 3.15.  Alternative Approaches
##

# p. 141
# example 3.6

data(d.sp8099)

d.lsp <- log(1+d.sp8099)

fit1.dlsp <- arima(d.lsp)
fit2.dlsp <- arima(d.lsp, c(0,0,1))

library(fGarch)
#fit1.dlsp <- garchFit(~1, data=d.lsp)
#fit2.dlsp <- garchFit(~arma(0, 1), data=d.lsp)

fit3.dlsp <- garchFit(~garch(1,1), data=d.lsp)
fit3.dlsp
# rather different from the book ...

# p. 142
# Figure 3.11

plot(fit3.dlsp@sigma.t) # ?????



# p. 143
# Figure 3.12
#  data not available


