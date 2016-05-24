asypow.sig <- function(asypow.obj, sample.size, power) {
#----------------------------------------------------------------------
#
#          Calculates the significance level of the test.
#
# asypow.obj : The object returned from asypow.noncent
#
# sample.size : The sample size of the test
#
# power : The power of the test.
#
#
# RETURNS the significance level
#
#----------------------------------------------------------------------

      n <- length(sample.size)
      m <- length(power)
      if (n !=1 & m!= 1 & n != m)
                stop("lengths of sample.size and power must match")
      if (n == 1) sample.size <- rep(sample.size, length=m)
      if (m == 1) power <- rep(power, length=n)

      w <- asypow.obj$w * sample.size

      crit <- rep(0,n)

      for(i in 1:n) {
	crit[i] <- cdfchn(2,p=1-power[i], df=asypow.obj$df, pnonc=w[i])$x
      }

      return(1-pchisq(crit,1))
}
