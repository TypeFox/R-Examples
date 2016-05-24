#############################################################################
# inverse of gamma  CDF  - inverse of F(x|alpha, delta) with respect to x   #
#############################################################################

inxncgamma   <-   function(p, alpha, delta)
{
   errortol <- 1e-12
   maxitr <- 50000
   d      <- delta
   k      <- ceiling(d / 2)
   a      <- alpha + k
   z      <- qnorm(p)  
   x0     <- ((a+4*d)*(z+((a+2*d)^2/(a+4*d)-1)^0.5)^2)/(a+2*d)
   xn     <- x0
   it     <-  1
   convnewton <- FALSE
   while (!convnewton)
   {
     x      <- xn
     gamac  <- pgamma(x,a)
     gamad  <- gamac
     gxd    <- exp(a*log(x)-x-lgamma(a+1))
     gxc    <- gxd * a / x
     hxc    <- gxd * a / x
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
         hxc    <- hxc * x / (a + i - 1)
         ppoic  <- ppoic * (d / 2) / (k + i)
         cdf    <- cdf + ppoic * gamac
         g      <- g + ppoic * hxc
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
            hxd    <- hxd * (a - i) / x
            ppoid  <- ppoid * (k - i + 1) / (d / 2)
            cdf    <- cdf + ppoid * gamad
            g      <- g + ppoid * hxd      
            remain <- remain - ppoid
            if ((remain <= errortol) | (i > maxitr)) convergiu <- TRUE
            i      <- i + 1
         }   
      } # achieved convergence
      if (x - (cdf - p) / g <= 0)
      {
        xn   <- x / 2
      }else
      {
        xn   <- x - (cdf - p) / g
      }
      if ((abs(xn - x) <= x * errortol) | (it > maxitr)) convnewton <- TRUE
      it <- it + 1
   } # end while convnewton
   return(xn)
} 
