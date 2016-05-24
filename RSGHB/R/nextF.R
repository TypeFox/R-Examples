nextF <- function(p,f,b,env)
{
     
     fnew <- f + sqrt(env$rhoF) * rnorm(env$gFIV,0,1)
     
     pnew <- env$likelihood(fnew,b,env)
     
     #r    <- prod(pnew/p)
     
     r    <- sum(log(pnew)-log(p))
          
     # ind  <- ( runif(1) <= r )
     ind  <- ( log(runif(1)) <= r )
          
     pnew <- pnew * ind + p * ( 1 - ind )
     
     fnew <- as.vector(t(fnew %*% t(ind)) + t(f %*% t( 1 - ind )))

     return(list(fnew,pnew))
}
