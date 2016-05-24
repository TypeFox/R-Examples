"cRR" <-
function (x)
 {
 rb = rbar(x)
 n <- length (x$u[!(is.na(x$u))])
 u <- x$u
 if (n == 0)
   { c <- 1 
     vc <- 0}
 else { c <- sqrt((1-u^2)*rb^2+u^2) 
 vc <- var(c, na.rm=TRUE)*(n-1)/n }
 mc <- mean(c, na.rm=TRUE)
 out <- list(mc, vc)
 return(out)
   }

