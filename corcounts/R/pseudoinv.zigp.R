`pseudoinv.zigp` <-
function(u, mu.x, phi.x=1, omega.x=0) {
     n <- length(u)
     x <- double(n)
     upper <- min(max(u),0.9999999)
     s <- double(1000)
     #P(X=0)
     p <- omega.x + (1-omega.x) * exp(-mu.x/phi.x)
     s[1] <- p
     if (upper > 0) {
        rekursive <- FALSE
        i <- 1
        while (s[i] < upper) {
          #P(X=x)
          if (rekursive==FALSE) {
            p <- (1-omega.x)*mu.x*(mu.x+(phi.x-1)*i)^(i-1)/exp(lgamma(i+1))*
                 phi.x^(-i)*exp(-1/phi.x*(mu.x+(phi.x-1)*i)) }
          if (p==Inf) {
            rekursive <- TRUE
            log.p.alt <- log( (1-omega.x)*mu.x*(mu.x+(phi.x-1)*(i-1))^(i-2)/exp(lgamma(i-1+1))*
                         phi.x^(-(i-1))*exp(-1/phi.x*(mu.x+(phi.x-1)*(i-1))) )
            }
          if (rekursive==TRUE) {
            log.p <- log( (mu.x+(i-1)*(phi.x-1))/(phi.x*i)*
                     (1+(phi.x-1)/(mu.x+(i-1)*(phi.x-1)))^(i-1)*
                     exp(1/phi.x-1) ) + log.p.alt
            log.p.alt <- log.p
            p <- exp(log.p)
          }
          if (ceiling(i/50)==floor(i/50)) {
            temp <- double(50)
            s <- c(s,temp)
          }
          s[i+1] <- s[i] + p
          i <- i+1
        }
     }

    for (j in 1:n) {
       i <- 1
       while (min(u[j],0.9999999) > s[i]) {
          i <- i+1
       }
       x[j] <- i-1
     }
   return(x)

}

