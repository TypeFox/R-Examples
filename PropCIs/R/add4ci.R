add4ci <-
function(x,n,conf.level){
   ptilde = (x+2)/(n+4)
   z = abs(qnorm((1-conf.level)/2))
   stderr = sqrt(ptilde * (1-ptilde)/(n+4))
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

