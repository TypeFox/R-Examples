betaval <-function(mn,sdev,fx=runif(1))
{
   if(mn>1 || mn<0){stop("Please select a mean beta between 0 and 1")} 
   # Original MATLAB script by Morris & Doak (2002: 277- 278),
   # adapted to R by Patrick Nantel, 20 June 2005.
   if(sdev == 0){bb <- mn}
   else
   {
      toler <- 0.0001 # this is tolerance of answer: how close
      # the CDF value of the answer must be to the input value (Fx)

      # this checks that the input mean and st. deviation
      # are possible for a beta.
      vari <- sdev^2
      if (vari >= (1 - mn) * mn)
      {
         stop('Standard deviation too high for beta distribution')
      }
      # start with a beginning guess x; the use of runif
      # adds wiggle to the search start to avoid pathologies
      vv <- mn * ((mn * (1 - mn)/(vari)) - 1) # calculate the beta parameters
      ww <- (1 - mn) * ((mn * (1 - mn)/(vari)) - 1)
      upval <- 1; lowval <- 0; x <- 0.5 + 0.02 * runif(1)
      i <- pbeta(x,vv,ww) # find the CDF value for x
      # the following while loop searches for ever better
      # values of x, until the value has a CDF within the
      # toler of Fx (unless the value
      # is very close to 0 or 1, which will also terminate
      # the search)
      while((toler < abs(i - fx)) & (x > 1e-6) & ((1 - x) > 1e-6))
      {
         if (fx > i)
         {
            lowval <- x
            x <- (upval + lowval)/2
         }
         else
         {
             upval <- x
             x <- (upval + lowval)/2
         } 
         i <- pbeta(x,vv,ww)
      } 
      # This makes values of x somewhat random to eliminate
      # pathologies when variance is very small or large.
      # It also truncates values of x, with the
      # smallest values equal to toler and the biggest
      # equal to 1 - toler.
      bbb <- x + toler * 0.1 * (0.5 - runif(1))
      if (bbb < toler)    bbb <- toler
      if (bbb > 1)        bbb <- 1 - toler
      bb <- bbb
   }   
   bb
}

