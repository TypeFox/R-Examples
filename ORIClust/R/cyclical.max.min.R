cyclical.max.min <-
function(data,x,n.rep,max1,min1){ 
   k <- length(n.rep)
   data <- as.numeric(data)
   n <- sum(n.rep)
   xa <- x
   waa <- n.rep
   xaa <- xa*waa 
   xpb <- rep(0,k)

   h1 <- isoincre(xaa[1:(max1-1)],waa[1:(max1-1)])

   if((max1+1)<=(min1-1))
      h2 <- isodecre(xaa[(max1+1):(min1-1)],waa[(max1+1):(min1-1)])
   if((max1+1)>(min1-1))
      h2 <- NULL


   h3 <- isoincre(xaa[(min1+1):k],waa[(min1+1):k])

   peak1 <- max(c(xa[max1],h1[length(h1)],h2[1]))
   peak2 <- min(c(xa[min1],h2[length(h2)],h3[1]))

   xpb <- c(h1,peak1,h2,peak2,h3)

   is <- xpb
 
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

