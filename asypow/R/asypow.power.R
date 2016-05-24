asypow.power <- function(asypow.obj, sample.size, significance) {
#----------------------------------------------------------------------
#
#               Calculates the power of a test
#
# asypow.obj : The object returned from asypow.noncent
#
# sample.size : The sample size of the study.
#
# significance: The significance level of the test.
#
#
# RETURNS the power of the test.
#
#----------------------------------------------------------------------

      n <- length(significance)
      m <- length(sample.size)
      if (n !=1 & m!= 1 & n != m)
                stop("lengths of significance and sample.size must match")
      if (n == 1) significance <- rep(significance,length=m)
      if (m == 1) sample.size <- rep(sample.size,length=n)

      crit <- qchisq(1-significance, asypow.obj$df)
      w <- asypow.obj$w * sample.size

      power <- rep(0,n)

      for(i in 1:n) {
	power[i] <- 1-cdfchn(1,x=crit[i],df=asypow.obj$df,pnonc=w[i])$p
      }
      return(power)
}

