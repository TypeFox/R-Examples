qda <-
function(x,y) {
     K <- max(y)
     p <- ncol(x)
     n <- nrow(x)
     m <- NULL  
     v <- array(0,c(K,p,p))
     ck <- NULL 
     phat <- table(y)/n  
     for(k in 1:K) {
          xk <- x[y==k,]
          m <- rbind(m,apply(xk,2,mean))  
          nk <- nrow(xk)
          v[k,,] <- var(xk)*(nk-1)/nk
          ck <- c(ck,-log(det(v[k,,]))+2*log(phat[k]))
     }

     list(Mean = m,Sigma = v, c = ck)
}
