#############################################################################
# noncentral pdf        -  density function                                 #
# f(x|alpha, delta)                                                         #
#############################################################################

dncgamma   <-   function(x, alpha, delta)
{
   errortol <- 1e-12
   maxitr <- 50000
   if (delta == 0)
   {
      if (x < 0) fdp <- 0 else 
      if (x==0)
      {
         if (alpha == 1) fdp <- 1 else
         if (alpha < 1) fdp <- Inf else
         fdp <- 0         
      } else
      fdp <- exp((alpha-1) * log(x) - x - lgamma(alpha)) 
      return(fdp) 
   } else
   {
     d      <- delta
     k      <- ceiling(d / 2)
     a      <- alpha + k
     hxc    <- 1/gamma(a+1) * exp(-x) * x^a
     if (x == 0) hxc <- 0 else
     hxc    <- hxc * a / x
     hxd    <- hxc
     ppoic  <- dpois(k, d / 2)
     ppoid  <- ppoic
     remain <- 1 - ppoic
     g      <- ppoic * hxc
     i      <- 1
     convergiu <- FALSE
     while (!convergiu)
     {
        hxc    <- hxc * x/(a + i -1)
        ppoic  <- ppoic * (d/2)/(k + i)
        g      <- g + ppoic * hxc
        error  <- remain * g
        remain <- remain - ppoic
        if (i > k)
        {
          if ((error <= errortol) | (i > maxitr)) convergiu <- TRUE
          i   <- i + 1
        } else 
        {
           hxd    <- hxd * (a - i) / x
           ppoid  <- ppoid * (k - i + 1) / (d / 2)
           g      <- g + ppoid * hxd      
           remain <- remain - ppoid
           if ((remain <= errortol) | (i > maxitr)) convergiu <- TRUE
           i      <- i + 1
         }   
     } # achieved convergence
     return(g)
   }
}
