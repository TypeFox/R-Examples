ccsize <- function(n,q,pD,p1,alpha,theta,power=NULL,verbose=FALSE)
{
   if (p1<0 | p1>1) stop("p1 should lie in [0,1]")
   if (pD<0 | pD>1) stop("pD should lie in [0,1]")
   if (q<0 | q>1) stop("q is should lie in [0,1]")
   p2 <- 1 - p1
   if(is.null(power))
# equation (5)
   {
     z_alpha <- qnorm(alpha)
     z <- z_alpha + sqrt(n) * theta * sqrt(p1 * p2 / (1 / pD + (1 / q - 1)))
     invisible(pnorm(z))
   }
   else
# equation (6)
   {
     nb <- -999
     z_alpha <- qnorm(alpha, lower.tail=FALSE)
     z_beta <- qnorm(power)
     theta_lon <- (z_alpha + z_beta) / sqrt(p1 * p2 * pD)
     d <- (theta / theta_lon)^2 - (1 - pD) / n
     if (d <= 0 & verbose) cat("bad hazard ratio =", exp(theta), "\n")
     else 
     {
       nb <- ceiling(pD / d)
       if (nb > n)
       {
         nb <- -999
         if (verbose) cat("bad subcohort size", nb, "\n")
       }
     }
     invisible(nb)
   }
}
