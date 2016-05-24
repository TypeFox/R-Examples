moment<-function(clz, vbz, output = "phi"){
     r <- length(clz)
     soma <- sum(vbz)
     mid <- clz
     mid[1] <- clz[1]-((clz[2]-clz[1])/2)
     for (i in 2:r) {mid[i] <- clz[i]-((clz[i]-clz[i-1])/2)}
     
     if(output == "metric"){
       m<-((1/2^mid) * 1000)
       m1<-log(m)
       p2 <- vbz * m1
       mean<-exp(sum(p2)/soma)
       std<-sum(vbz*((m1-log(mean))^2))
       sort<-exp(sqrt(std/soma))
       std.skew<-sum(vbz*((m1-log(mean))^3))
       skew<-std.skew / (soma* (log (sort)^3) )
       std.kurt<-sum(vbz*((m1-log(mean))^4))
       kurt<-std.kurt / (soma* (log (sort)^4) )
       std.A5<-sum(vbz*((m1-log(mean))^5))
       A5<-std.A5 / (soma* (log (sort)^5) )
       std.A6<-sum(vbz*((m1-log(mean))^6))
       A6<-std.A6 / (soma* (log (sort)^6) )
       }
     if (output=="phi"){
       vfz <- vbz*100/soma
       pz <- vfz
       pz[1] <- vfz[1]
       for (i in 2:r) pz[i] <- pz[i-1] + vfz[i]
       pz <- mid * vbz
       mean <- sum(pz)/soma
       std<-mid-mean
       pz <- vbz * std^2
       sort1 <- sum(pz)/(soma)
       sort <- sqrt(sum(pz)/(soma))
       pz <- vbz * std^3
       skew1 <- sum(pz)/(soma)
       skew <- skew1/(exp(1.5*log(sort1)))
       pz <- vbz * std^4
       kurt1 <- sum(pz)/(soma)
       kurt <- kurt1/(exp(2*log(sort1)))
       pz <- vbz * std^5
       A5 <- sum(pz)/(soma)
       A5 <- A5/(exp(2.5*log(sort1)))
       pz <- vbz * std^6
       A6 <- sum(pz)/(soma)
       A6 <- A6/(exp(3*log(sort1)))
       }
  table<-list(mean=mean, sort=sort, skew=skew, kurt=kurt, A5=A5, A6=A6)
  return(table)
   }
