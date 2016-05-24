"F0w.s" <-
function(u) {
# Fo+(u) (the function adaptw uses a Fortran version)
if (u > 20) return(1)
if (u <=1) z <- 0 else
 {if (u > 1.5) Tl <- uniroot(rhow,lower=-u,upper=-u+1.5,const=u)$root else
               Tl <- uniroot(rhow,lower=-u,upper=0,const=u)$root 
               Tu <- uniroot(rhow,lower=log(u),upper=u,const=u)$root
z <- pweibull(exp(Tu),shape=1)-pweibull(exp(Tl),shape=1)}; z}

