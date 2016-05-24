mplik.wb.s <-
function(par,Y,X,delta){
         y <- log(Y) 
     sigma <- exp(par[1])
       phi <- matrix(par[-1],ncol=1)
         p <- nrow(phi)         
         r <- sum(delta)
         z <- as.vector((y-X%*%phi)/sigma)
         Z <- diag(exp(z))
         V <- matrix(NA,nrow=length(y),ncol=p)
   for(i in 1:length(y)){
           if(delta[i]==1){
              V[i,] <- X[i,]
          }else{
              V[i,] <- 0   
          }                  
        }
loglik = -r*log(sigma)+sum((y-X%*%phi)*delta)/sigma-sum(exp((y-X%*%phi)/sigma))
 fn = loglik+ p*log(sigma)+.5*log(abs(det(t(X)%*%Z%*%X)))-log(abs(det(t(X)%*%Z%*%V)))
return(-fn)
}
