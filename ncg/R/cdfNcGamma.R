#############################################################################
# noncentral gamma cdf                                                      #
# F(x|alpha, delta)                                                         #
#############################################################################

pncgamma   <-   function(x, alpha, delta)
{
   errortol <- 1e-12
   if (x < 0) x <- 0  
   maxitr <- 50000
   d      <- delta
   k      <- ceiling(d / 2)
   a      <- alpha + k
   gamac  <- pgamma(x, a)
   gamad  <- gamac
   gxd    <- exp(a * log(x) - x - lgamma(a+1))
   if (x == 0) gxc <- 0 else 
     gxc    <- gxd * a / x
   ppoic  <- dpois(k, d / 2)
   ppoid  <- ppoic
   remain <- 1 - ppoic
   cdf    <- ppoic * gamac
   i       <-  1
   convergiu <- FALSE
   while (!convergiu)
   {
      gxc    <- gxc * x / (a + i - 1)
      gamac  <- gamac - gxc
      ppoic  <- ppoic * (d / 2) / (k + i)
      cdf    <- cdf + ppoic * gamac
      error  <- remain * gamac
      remain <- remain - ppoic
      if (i > k)
      {
        if ((error <= errortol) | (i > maxitr)) convergiu <- TRUE
         i   <- i + 1
      } else 
      {
        gxd    <- gxd * (a + 1 - i) / x
        gamad  <- gamad + gxd
        ppoid  <- ppoid * (k - i + 1) / (d / 2)
        cdf    <- cdf + ppoid * gamad      
        remain <- remain - ppoid
        if ((remain <= errortol) | (i > maxitr)) convergiu <- TRUE
        i      <- i + 1
      }   
   } # achieved convergence 
   return(cdf)   
}
