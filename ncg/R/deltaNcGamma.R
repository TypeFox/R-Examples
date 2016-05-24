indncgamma<- function(x, alpha, p)
{
   errortol <- 1e-12
   maxitr <- 50000
   dn     <- d0(x,p,alpha)
   it     <-  1
   convnewton <- FALSE
   while (!convnewton)
   {
      d      <- dn      
      k      <- ceiling(d / 2)
      a      <- alpha + k
      gamac  <- pgamma(x,a)
      gamad  <- gamac  
      gxd    <- exp(a*log(x)-x-lgamma(a+1))
      gxc    <- gxd * a / x
      hxc    <- gxd
      hxd    <- hxc
      ppoic  <- dpois(k, d / 2)
      ppoid  <- ppoic
      remain <- 1 - ppoic
      cdf    <- ppoic * gamac
      g      <- ppoic * hxc
      i      <-  1
      convergiu <- FALSE
      while (!convergiu)
      {
         gxc    <- gxc * x / (a + i - 1)
         gamac  <- gamac - gxc
         hxc    <- hxc * x / (a + i)
         ppoic  <- ppoic * (d / 2) / (k + i)
         cdf    <- cdf + ppoic * gamac
         g      <- g + ppoic * hxc
         error  <- remain * gamac
         remain <- remain - ppoic
         if (i > k)
         {
            if ((error <= errortol) | (i > maxitr)) convergiu <- TRUE
            i   <- i + 1
         }else 
         {
            gxd    <- gxd * (a + 1 - i) / x
            gamad  <- gamad + gxd
            hxd    <- hxd * (a + 1 - i) / x
            ppoid  <- ppoid * (k - i + 1) / (d / 2)
            cdf    <- cdf + ppoid * gamad
            g      <- g + ppoid * hxd      
            remain <- remain - ppoid
            if ((remain <= errortol) | (i > maxitr)) convergiu <- TRUE
            i      <- i + 1
         }   
      } # achieved convergence
      if (d + 2*(cdf - p) / g <= 0)
      {
         dn   <- d / 2
      } else
      {
         dn   <- d + 2*(cdf - p) / g
      }
      if ((abs(dn - d) <= d * errortol) | (it > maxitr)) convnewton <- TRUE
      it <- it + 1
   } # end while convnewton
   return(dn)
}
