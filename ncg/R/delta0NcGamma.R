#############################################################################
# inverse of gamma CDF - inverse of F(x|alpha, delta) with respect to delta #
#############################################################################

###initial value of the noncentrality parameter -  d0###
d0<- function(x, p, alpha)
{
   dn0   <- 0.5
   dn1   <- 1
   errortol<- 1e-12
   nmax  <- 1000
   z     <- qnorm(p)
   a     <- alpha
   fdn0  <- ((a+4*dn0)*(z+((a+2*dn0)^2/(a+4*dn0)-1)^0.5)^2)/(a+2*dn0)-x
   fdn1  <- ((a+4*dn1)*(z+((a+2*dn1)^2/(a+4*dn1)-1)^0.5)^2)/(a+2*dn1)-x
   it    <- 1
   convd <- FALSE
   while(!convd)
   {
      dn2   <- dn1 - fdn1 * (dn1 - dn0) / (fdn1 - fdn0)
      if (dn2=="NaN") break
      if (dn2 < 0)   dn2 <- dn1 / 2
      dn0   <- dn1
      fdn0  <- fdn1
      dn1   <- dn2
      fdn1  <- ((a+4*dn1)*(z+((a+2*dn1)^2/(a+4*dn1)-1)^0.5)^2)/(a+2*dn1)-x
      if ((abs(dn2 - dn0) <= errortol * dn0) | (it > nmax)) convd <- TRUE
      it <- it + 1
   } #end while convd
   if (dn2=="NaN") return(1) else return(dn2) # protection for poor result
}
