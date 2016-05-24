### ch. 4.  Nonlinear models and their applications 
###
### 
### Ruey S. Tsay (2005)
### Analysis of Financial Time Series, 2nd ed.
### (Wiley)
###
### 
library(FinTS)
# p. 154

# p. 156 
##
## sec. 4.1.  Nonlinear models  
##

# sec. 4.1.1.  Bilinear model 

# p. 157
# Example 4.1.
data(m.ibmvwewsp2603)
str(m.ibmvwewsp2603)

ew2697 <- window(m.ibmvwewsp2603[, "VW"], end=yearmon(1997+11/12))
str(ew2697)
# correct number of observations.

#??? Need code for bilinear estimation.

# arch(3)
library(tseries)

# Best I currently know in R (2008.03.09) 
# (1) Fit an AR(3) 
# (2) Fit an ARCH(3) to the residuals:

(ew.ar3 <- arima(as.numeric(ew2697), c(3, 0, 0), fixed=c(NA, 0, NA, NA)))
(ew.ar3.arch3 <- garch(resid(ew.ar3), order=c(0, 3)))

# poor match with the numbers on p. 157

# sec. 4.1.2.  Threshold autoregression (TAR)

# p. 158
# Figure 4.1.  Simulated 2-regime TAR(1) Series

# p. 159
# Example 4.2.  US monthly unemployment
data(m.unrate)
str(m.unrate)

(unrateARIMA <- arima(m.unrate, c(2, 1, 2),
                     seasonal=list(order=c(1,0,1), period=12)))
# good match to (4.10) 
str(unrateARIMA)
sqrt(unrateARIMA$sigma2)

# p. 160
plot(m.unrate, xlab='year')
plot(m.unrate, type='b', xlab='year')

AutocorTest(resid(unrateARIMA), 12)
AutocorTest(resid(unrateARIMA), 24)
# close but not exactly the numbers in the book.

# Fit in the book estimes AR lags
# 2, 3, 4, 12 with y[t-1]<=0.1
# 2, 3,    12 with y[t-1]> 0.1

# 'setar' in tsDyn will estimate TAR models
# but will not allow parameters to be fixed like this

library(tsDyn)

quantile(as.numeric(diff(m.unrate)))
#  0%  25%  50%  75% 100% 
#-1.5 -0.1  0.0  0.1  1.3

data(m.unrate) 
unTAR <- setar(diff(m.unrate), 12, th=0.1)
unTAR
summary(unTAR)

# p. 161
# example 4.3
data(d.ibmvwewsp6203)
str(d.ibmvwewsp6203)
# correct number of observations

quantile(coredata(d.ibmvwewsp6203[, "IBM"]))
#         0%        25%        50%        75%       100% 
#-0.2296300 -0.0083800  0.0000000  0.0088075  0.1316400 
l.ibm <- log(1+d.ibmvwewsp6203[, "IBM"])
quantile(coredata(l.ibm))
#         0%         25%         50%         75%        100% 
#-0.26088436 -0.00841531  0.00000000  0.00876894  0.12366791

library(fGarch) 
fit2.garch1.1 <- garchFit(~arma(2,0)+garch(1,1), coredata(l.ibm) ) 
fit2.garch1.1
                          
