scoreci.mp <-
function(b,c,n,conf.level)
{
   pa = 2*n
   z = qnorm(1-(1-conf.level)/2)

   if(c == n) {ul = 1}
   else{
     proot = (c-b)/n
     dp = 1-proot
     niter = 1
     while(niter <= 50){
       dp = 0.5*dp
       up2 = proot+dp
       pb = - b - c + (2*n-c+b)*up2
       pc = -b*up2*(1-up2)
       q21 = (sqrt(pb^2-4*pa*pc)-pb)/(2*pa)
       score = (c-b-n*up2)/sqrt(n*(2*q21+up2*(1-up2)))
       if(abs(score)<z){ proot = up2 }
       niter=niter+1
       if((dp<0.0000001) || (abs(z-score)<.000001)){
       niter=51
       ul=up2
       }
      }
   }

   if(b == n) {ll = -1}
   else{
     proot = (c-b)/n
     dp = 1+proot
     niter = 1
     while(niter <= 50){
       dp = 0.5*dp
       low2 = proot-dp
       pb = - b - c + (2*n-c+b)*low2
       pc = -b*low2*(1-low2)
       q21 = (sqrt(pb^2-4*pa*pc)-pb)/(2*pa)
       score = (c-b-n*low2)/sqrt(n*(2*q21+low2*(1-low2)))
       if(abs(score) < z){proot = low2}
       niter = niter+1
       if((dp<0.0000001) || (abs(z-score)<.000001)){
       ll = low2
       niter = 51
       }
      }
  }
   cint <- c(ll, ul)
   attr(cint, "conf.level") <- conf.level
   rval <- list(conf.int = cint)
   class(rval) <- "htest"
   return(rval)
}

