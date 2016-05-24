ab <- function(type="power",n=25000,a=0.15,sa=0.01,b=log(1.19),sb=0.01,alpha=0.05,fold=1)
{
   ab <- a*b
   s <- sqrt(a^2*sb^2+b^2*sa^2)
   z <- ab/s
   if (type=="power") 
   {
      x2 <- z^2
      x <- qchisq(alpha,1,lower.tail=FALSE)
      power <- pchisq(x,1,ncp=x2*fold,lower.tail=FALSE)
      cat(fold*n, ",", power, "\n")
      stats <- c(fold*n,power)
   } else if (type=="test") stats <- c(z,2*pnorm(-abs(z)))
   else stop("Invalid option")
   invisible(stats)
}

# 10-11-2009 Modified from Stata code
