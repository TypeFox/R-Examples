### ch. 11 State Space Models and Kalman Filter
###
### 
### Ruey S. Tsay (2005)
### Analysis of Financial Time Series, 2nd ed.
### (Wiley)
###
###
library(FinTS)

#  p. 490
##
## Sec. 11.1 Local Trend Model
## 

#  p. 492
##
## Example 11.1 
##

data(aa.3rv)

(arimaest <- ARIMA(as.numeric(log(aa.3rv[,"X10m"])),
                   order = c(0,1,1), Box.test.lag=12))
names(arimaest)
sqrt(arimaest$sigma2)
arimaest$Box.test

#=================================================
# ESTMATION ABOVE GIVES S.E OF THETA AS .0397. 
# RUEY HAS .029 ON BOTTOM OF PG 492
#+++++++++++++++++++++++++++++++++++++++++++++++++

(boxtestressq<-Box.test((arimaest$residuals)^2,lag=12,type="Ljung-Box"))

#=======================================================================
# SG POINTED OUT THAT DIFFERENCE IN FIRST PVALUE IS DUE TO USING 12
# FOR DF RATHER THAN 11. 

#        Box-Ljung test
#
#data:  arimaest$residuals 
#X-squared = 12.4128, df = 12, p-value = 0.4131
#
#
#        Box-Ljung test
#data:  (arimaest$residuals)^2 
#X-squared = 8.1606, df = 12, p-value = 0.7725

# 1-pchisq(8.16,12)
#[1] 0.7725052

# 1-pchisq(12.41,12)
#[1] 0.4133387
#==========================================================================

# p. 493

#===================================================================
# SG:SIGN OF MA COEFF CONVENTION IS OPPOSITE IN R
# I USE AS.NUMERIC TO TAKE NAME OFF BUT MAYBE THERE IS  A BETTER WAY
#====================================================================

(thetacoeff<- as.numeric(-1.0*arimaest$coef))
(sigmasquareda <- arimaest$sigma2)

(sigmaepsilon<- sqrt(sigmasquareda*thetacoeff))
(sigmaeta <- sqrt((1+thetacoeff^2)*sigmasquareda - 2*sigmaepsilon^2))

# Figure 11.1. Time plot of the logarithms of intraday realized 
#  	       volatility of Alcoa stock from January 2, 2003 to May 7, 2004. 
#	       The realized volatility is computed from the intraday 10-minute 
#	       log returns measured in percentage.

plot(as.numeric(log(aa.3rv[,"X10m"])), type="l",xlab="day", ylab="rv")

# Sec. 11.1.2  Kalman Filter

# p. 496
##
## Example 11.1 Continued
##

arima(diff(as.numeric(log(aa.3rv[,"X10m"]))),order = c(0,0,1),include.mean= FALSE)

library(dlm)

#======
# NEW
buildRwpn <- function(x) {
paramlist<- dlmModPoly(1, dV = exp(x[1]), dW = exp(x[2]))
return(paramlist)
}
#=====

#=======
# NEW
(RwpnMLE <- dlmMLE(as.numeric(log(aa.3rv[,"X10m"])), rep(0,2), buildRwpn))
#=======

#=========
# MODIFIED-WRONG EARLIER
(alcoaMod <- dlmModPoly(1, dV = (exp(RwpnMLE$par[1])), dW = (exp(RwpnMLE$par[2]))))
#==========

(alcoaFilt <- dlmFilter(as.numeric(log(aa.3rv[,"X10m"])),alcoaMod))
(alcoaSmooth <- dlmSmooth(alcoaFilt))

(alcoaCt <- with(alcoaFilt,dlmSvd2var(U.C,D.C)))
(alcoaZt <- with(alcoaSmooth,dlmSvd2var(U.S,D.S)))

#=========
#NEW
(resids <- residuals(alcoaFilt, type="raw")$res)
(stdresids <- residuals(alcoaFilt, type="standardized")$res)
#=========

#==============================================================
# NEW
# BOX TEST ON RESIDUALS AND RESIDUALS SQUARED BELOW GIVES DIFFERENT RESULT THAN BOOK
# NEED TO FIND OUT WHAT MOUT IN SSFPACK REPRESENTS
(boxteststdres<-Box.test(stdresids,lag=25,type="Ljung-Box"))
(boxteststdressq<-Box.test(stdresids^2,lag=25,type="Ljung-Box"))
#==============================================================

# p. 497
# Figure 11.2. Time plots of output of the Kalman Filter applied to the
#  	       daily realized log volatility of Alcoa stock based on the
#	       local trend state-space model: (a) the filtered state u_t|t
#	       and (b) the one-step ahead forecast error v_t.

op <- par(mfrow=c(2,1))
plot(alcoaFilt$m[-1],type="l",xlab="day",ylab="filtered state",main="(a) Filtered state variable")
# MODIFIED
plot(resids,type="l",xlab="day",ylab="v(t)",main="(b) Prediction error")
par(op)

# Sec. 11.1.4  State Smoothing

# p. 501
##
## Example 11.1 Continued
##

(mupper <- alcoaFilt$m + 1.96*sqrt(as.numeric(alcoaCt)))
(mlower <- alcoaFilt$m - 1.96*sqrt(as.numeric(alcoaCt)))
(supper <- alcoaSmooth$s + 1.96*sqrt(as.numeric(alcoaZt)))
(slower <- alcoaSmooth$s - 1.96*sqrt(as.numeric(alcoaZt)))

# p. 502
# Figure 11.3. Filtered state variabe u_t|t and its 95% confidence interval
#  	       for the daily log realized volatility of Alcoa stock returns based 
#	       on the fitted local trend state-space model.
	       
plot(alcoaFilt$m[-1],type="l",xlab="day",ylab="value",ylim=c(min(mlower[-1]- 0.1),2.5),lty=1)
lines(mupper[-1],type="l",xlab="day",ylab="value",lty=2)
lines(mlower[-1],type="l",xlab="day",ylab="value",lty=2)

# p. 502
# Figure 11.4. Smoothed state variable u_t|T and its 95% confidence interval 
#  	       for the daily log realized volatility of Alcoa stock returns
#	       based on the fitted local trend state-space model.
	       
plot(alcoaSmooth$s[-1],type="l",xlab="day",ylab="value",ylim=c(min(slower[-1]- 0.1),2.5),lty=1)
lines(supper[-1],type="l",xlab="day",ylab="value",lty=2)
lines(slower[-1],type="l",xlab="day",ylab="value",lty=2)


