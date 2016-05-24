up.down <-
function(data,x,n.rep,h){
   k <- length(n.rep)
   data <- as.numeric(data)
   n <- sum(n.rep)
   xa <- x
   waa <- n.rep
   xaa <- xa*waa
   xpb <- rep(0,k)
   em1 <- matrix(0,nrow=h,ncol=(k-h+1))
   for(s in 1:h)
      for(t in 1:(k-h+1))
         em1[s,t]<- sum(xaa[s:(t+h-1)])/sum(waa[s:(t+h-1)])
   
   xpb[h] <- max(em1)
   h1 <- isoincre(xaa[1:(h-1)],waa[1:(h-1)])
   h1[h1>xpb[h]] <- xpb[h]
   xpb[1:(h-1)] <- h1
   h2 <- isodecre(xaa[(h+1):k],waa[(h+1):k])
   h2[h2>xpb[h]] <- xpb[h]
   xpb[(h+1):k] <- h2
   sig <- 0;
   len1 <- 1;
   len2 <- 0;
   for(i in 1:k){
      muu <- xpb[i]
      len2 <- len2+n.rep[i]
      x1 <- data[len1:len2]
      len1 <- len1+n.rep[i]
      sig <- sig+sum((x1-muu)^2)
   }
   sig <- sig/n;
   like <- -0.5*n*log(2*pi)-0.5*n*log(sig)-0.5*n

   list( logelr=like, mu=xpb)
}

