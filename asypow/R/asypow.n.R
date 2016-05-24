asypow.n <- function(asypow.obj, power, significance) {
#----------------------------------------------------------------------
#
#          Calculates the sample size needed to obtain the
#                      desired power for a test.
#
# asypow.obj : The object returned from asypow.noncent
#
# power : The power of the test.
#
# significance: The significance level of the test.
#
#
# RETURNS the needed sample size.
#
#----------------------------------------------------------------------

      n <- length(significance)
      m <- length(power)
      if (n !=1 & m!= 1 & n != m)
                stop("lengths of significance and power must match")
      if (n == 1) significance <- rep(significance,length=m)
      if (m == 1) power <- rep(power,length=n)

      w <- rep(0,n)

      crit <- qchisq(1-significance,asypow.obj$df)

      for(i in 1:n) {
	 w[i] <- cdfchn(4, p=1-power[i], x=crit[i], df=asypow.obj$df)$pnonc
      }

      sample.size <- w / asypow.obj$w

      return(sample.size)
}
