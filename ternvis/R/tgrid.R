tgrid <-
function(n){
   p1 <- matrix(0,nrow=0.5*(n+1)*(n+2))
   p2 <- matrix(0,nrow=0.5*(n+1)*(n+2))
   p3 <- matrix(0,nrow=0.5*(n+1)*(n+2))
   c <- 0  
   for (i in 0:n){
      for (j in 0:(n-i)){
            c <- c + 1
	    p1[c] <- i/n
	    p2[c] <- j/n
	    p3[c] <- (n-i-j)/n
      }
   }
   out <- cbind(p1,p2,p3)
   out
}
