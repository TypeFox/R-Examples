skew<-function(x){
                     
          sk <- function(xx){
             n <- length(xx)
             mn <- mean(xx)
             dif.x <- xx - mn
             m2 <- sum(dif.x^2)/n
             m3 <- sum(dif.x^3)/n
             m4 <- sum(dif.x^4)/n
             b1 <- m3/(m2^(3/2))
             g1 <- (b1 * sqrt(n * (n - 1)))/(n - 2)
             g1
           }
          
          if(ncol(x)==1 || is.null(dim(x)))
             return(sk(x))
          else
           return(apply(x,2,sk))              
            
}
