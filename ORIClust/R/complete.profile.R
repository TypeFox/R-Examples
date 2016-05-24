complete.profile <-
function(data,x,n.rep){
   k <- length(n.rep)
   data <- as.numeric(data)
   n <- sum(n.rep)
   is <- x
   sig <- 0;
   len1 <- 1;
   len2 <- 0;
   for(i in 1:k){
       muu <- is[i]
       len2 <- len2+n.rep[i]
       x1 <- data[len1:len2]
       len1 <- len1+n.rep[i]
       sig <- sig+sum((x1-muu)^2)
   }
   sig <- sig/n;
   like <- -0.5*n*log(2*pi)-0.5*n*log(sig)-0.5*n
   list( logelr=like, mu=is)
}

