###################################################
### chunk number 1:
###################################################

######################################################################
# Workhorse computing the OutbreakP statistic.
# Alarm statistic at end time n is returned.
#
# Author:
# Michael Hoehle <hoehle@stat.uni-muenchen.de>
#
# R port of the Java code by Marianne Frisen & Linus Schioler from
# the CASE project. See https://smisvn.smi.se/case/
#
# For a manual on how to use the method see also
# http://www.hgu.gu.se/item.aspx?id=16857
#
# Date:
# 25 May 2010
#
# Parameters:
#  x -- the series with the counts
#
# Returns:
#  value of the alarm statistic at the end of the series x.
######################################################################

calc.outbreakP.statistic <- function(x) {
  #Length of the monitored series
  n <- length(x)
  #Index problem when converting java arrays to R arrays
  x <- c(0,x)

  #Initialization (not all parts might be needed)
  leftl <- numeric(n+1);
  y <- numeric(n+1);
  yhat <- numeric(n+1);
  sumwy <- numeric(n+1);
  sumwys <- numeric(n+1);
  sumw <- numeric(n+1);
  w <- numeric(n+1);
  meanl <- numeric(n+1);

  xbar <- 0
  meanl[1] = -Inf
  leftl[1] = 0


  for (i in 1:n) {
    #Initialize
    yhat[i+1] <- x[i+1];
    sumwy[i+1] <- x[i+1];
    sumw[i+1] <- 1;
    meanl[i+1] <- x[i+1];
    leftl[i+1] <- i;
    #Calculate mean (this is a sequential formula to calculate mean(x[1:i]))
    xbar=xbar+(x[i+1]-xbar)/i

    #Create plateaus
    while (meanl[i+1] <= meanl[ (leftl[i+1] - 1) + 1]) {
      #merge sets
      sumwy[i+1] = sumwy[i+1] + sumwy[(leftl[i+1] - 1)+1];
      sumw[i+1] = sumw[i+1] + sumw[(leftl[i+1] - 1)+1];
      meanl[i+1] = sumwy[i+1] / sumw[i+1];
      leftl[i+1] = leftl[(leftl[i+1] - 1)+1];
    }

    #calculate yhat
    for (j in leftl[i+1]:i) {
      yhat[j+1] = meanl[i+1];
    }
  }

  #Compute the statistic in case of a Poisson distribution
  alarm.stat <- 1
  for (j in seq_len(n)) {
    #Ensure 0/0 = 1 so we don't get NaNs
    div <- ifelse(yhat[j+1]==0 & xbar==0, 1, yhat[j+1]/xbar)
    alarm.stat <- alarm.stat * (div)^x[j+1]
  }
  return(alarm.stat)

## The above might cause NaN's in case of large numbers.
##  logalarm <- 0
##  for (j in 1:n) {
##     #Eqn (5) in Frisen et al paper in log form. However: it is undefined
##     #if \hat{\mu}^D(t) == 0 (it is a division by zero).
##     #We fix 0/0 = 1
##     if (xbar != 0) {
##       if (yhat[j+1] != 0) { #if \hat{\mu}^{C1} == 0 then
##         logalarm = logalarm + x[j+1] * (log(yhat[j+1]) - log(xbar))
##       }
##     } else {
##       if (yhat[j+1] != 0) {
##         stop("Division by zero in Eqn (5) of Frisen paper!")
##       }
##     }
##  }
##  #Done, return the value
##  return(exp(logalarm))
}

######################################################################
# The detection function in S3 style
######################################################################

algo.outbreakP <- function(disProgObj, control = list(range = range, k=100, ret=c("cases","value"),maxUpperboundCases=1e5)) {
  #Set threshold to some fixed value, i.e. 100
  if(is.null(control[["k",exact=TRUE]]))
    control$k <- 100

  #Set largest observed value to try as upperbound when numerically searching
  #for NNBA in case ret = "cases"
  if(is.null(control[["maxUpperboundCases",exact=TRUE]]))
    control$maxUpperboundCases <- 1e5

  #Which value to return in upperbound?
  control$ret <- match.arg(control$ret, c("value","cases"))

  #Initialize the necessary vectors
  alarm <- matrix(data = 0, nrow = length(control$range), ncol = 1)
  upperbound <- matrix(data = 0, nrow = length(control$range), ncol = 1)
  observed <- disProgObj$observed

  #Store results
  count <- 1
  for(i in control$range) {
    statistic <- calc.outbreakP.statistic( observed[seq_len(i)] )
    # store the results in the right order
    alarm[count] <- statistic > control$k

    #Find NNBA or just return value of the test statistic (faster)
    if (control$ret == "cases") {
      #If length is 1 no alarm can be given unless k<1
      if (i<=1) {
        upperbound[count] <- ifelse(control$k>=1, NA, 0)
      } else {
        if (is.nan(statistic)) { #if no decent statistic was computed.
          upperbound[count] <- NA
        } else {
          #Go up or down
          delta <- ifelse(alarm[count], -1, 1)
          #Initialize
          observedi <- observed[i]
          foundNNBA <- FALSE
          #Loop with modified last observation until alarm is caused (dx=1)
          #or until NO alarm is caused anymore (dx=-1)
          while ( ((delta == -1 & observedi > 0) | (delta == 1 & observedi < control$maxUpperboundCases)) & (!foundNNBA)) {
            observedi <- observedi + delta
            newObserved <- c(observed[seq_len(i-1)],observedi)
            statistic <- calc.outbreakP.statistic( newObserved )
            if (is.nan(statistic)) { #statistic produced a numeric overflow.
              observedi <- control$maxUpperboundCases
            } else {
              foundNNBA <- (statistic > control$k) == ifelse(alarm[count],FALSE,TRUE)
            }
          }
          upperbound[count] <- ifelse( foundNNBA, observedi + ifelse(alarm[count],1,0), NA)
        }
      }
    } else {
      upperbound[count] <- statistic
    }

    #Advance time index
    count <- count + 1
  }

  #Add name and data name to control object.
  control$name <- paste("outbreakP(",control$k,")",sep="")
  control$data <- paste(deparse(substitute(disProgObj)))

  # return alarm and upperbound vectors
  result <- list(alarm = alarm, upperbound = upperbound, disProgObj=disProgObj, control=control)

  class(result) = "survRes" # for surveillance system result
  return(result)
}



