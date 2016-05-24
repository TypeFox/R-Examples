# Huber's M-estimator; choice q=pchisq(k+1,k) gives maximum breakdown point

huber.M<-function( data, q=pchisq(3,2) )
{
  n <- dim(data)[1]
  k <- dim(data)[2]
  c <- qchisq(q,k)
  c2 <- pchisq(c,k+2)+(c/k)*(1-q)
  bdp <- min(1/c,1-k/c)
  MAX_ITER <- 250; # Max number of iteration
  EPS <- 1.0e-6;   # Iteration accuracy

# STARTING VALUES 

  t <- apply(data,2,mean)
  C <- cov(data)
  R<-chol(solve(C))

  iter <- 0
  d1<-1
  d2<-1

  while (d1 > EPS & d2 > EPS ) {       
# SCATTER STEP
      u<-NULL
      data2 <- sweep(data,2,t,"-")   
      s <- diag(data2%*%(t(R)%*%R)%*%t(data2))
      u[s<=c] <- (1/c2)
      u[s>c] <- (s[s>c]/(c/c2))^{-1}
      C <- R%*%(t(data2)%*%sweep(data2,1,u,"*")/n)%*%t(R)
      R0 <- chol(solve(C))
      R <- R0%*%R           
      d1 <- max(apply(abs(R0-diag(k)),1,sum))
# LOCATION STEP
      v <- NULL
      s <- diag(data2%*%(t(R)%*%R)%*%t(data2)) 
      v[s<=c] <- 1
      v[s>c] <- sqrt(c)/sqrt(s[s>c])     
      h <- apply(sweep(data2,1,v,"*"),2,mean)/mean(v);
      t <- t + h;
      d2 <-sqrt(t(h-t)%*%(h-t))
      iter <- iter+1     
      if (iter > MAX_ITER) break 
  }           

  C <- t(R)%*%R
  C <- solve(C) 
list (loc=t,cov=C,iter=iter)
}   
