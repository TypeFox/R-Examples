lda <-
function(x,y) {
     if(is.vector(x)) {x <- matrix(x,ncol=1)}
     K <- max(y)
     p <- ncol(x)
     n <- nrow(x)
     m <- NULL  
     v <- matrix(0,ncol=p,nrow=p) 
     for(k in 1:K) {
          xk <- x[y==k,]
          if(is.vector(xk)) {xk <- matrix(xk,ncol=1)}
          m <- rbind(m,apply(xk,2,mean))  
          v <- v + var(xk)*(nrow(xk)-1)
     }
     v <- v/n
     phat <- table(y)/n  

     ck <- NULL  
     ak <- NULL  
     vi <- solve(v)  
     for(k in 1:K) {
          c0 <- -(1/2)*(m[k,]%*%vi%*%m[k,]-m[K,]%*%vi%*%m[K,])
                                                          +log(phat[k]/phat[K])
          ck <- c(ck,c0)
          a0 <- vi%*%(m[k,]-m[K,])
          ak <- cbind(ak,a0)
     }
     list(a = ak, c = ck)
}
