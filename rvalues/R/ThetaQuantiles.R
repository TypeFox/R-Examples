
ThetaQuantiles <- function(FF, alpha.grid, lbd, ubd) {
    
    ## FF should be a function produced from approxfun or stepfun
     
    if(is.stepfun(FF)) {
         kk <- knots(FF)
         kk <- rev(kk)
         nknots <- length(kk)
         
         Fvals <- c(1 - FF(kk),1)
         
         qtheta <- rep(0, length(alpha.grid))
         for(i in 1:nknots) {
             qtheta[(alpha.grid >= Fvals[i]) & (alpha.grid < Fvals[i+1])] <- kk[i]
         }
    }
    else {
         ff <- function(x, alpha)  {
              return(1 - FF(x) - alpha)
         }
         ngrid <- length(alpha.grid)
         qtheta <- mroot(ff, lower=rep(lbd, ngrid), upper=rep(ubd, ngrid), 
                         alpha = alpha.grid)$root 

         #for(i in 1:ngrid)  {
         #    qtheta[i] <- uniroot(ff, lower = lbd, upper = ubd, alpha = alpha.grid[i])$root
         #}
    }    
    return(qtheta)   
}
