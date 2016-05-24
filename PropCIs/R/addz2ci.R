addz2ci <-
function(x,n,conf.level){
   z = abs(qnorm((1-conf.level)/2))
   tr = z^2     #the number of trials added
   suc = tr/2   #the number of successes added
   ptilde = (x+suc)/(n+tr)
   stderr = sqrt(ptilde * (1-ptilde)/(n+tr))
   ul = ptilde + z * stderr
   ll = ptilde - z * stderr
   if(ll < 0) ll = 0
   if(ul > 1) ul = 1
   cint <- c(ll, ul)
   attr(cint, "conf.level") <- conf.level
   rval <- list(conf.int = cint, estimate = ptilde)
   class(rval) <- "htest"
   return(rval)
}

