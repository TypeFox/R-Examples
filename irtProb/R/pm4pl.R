`pm4pl` <-
function(theta=0,S=0,C=0,D=0,s=1/1.702,b=0,c=0,d=1,lower.tail=TRUE,log.p=FALSE) {
 result <- (C+c) + ((d-D) - (C+c))/(1+exp(-(theta-b)/sqrt(s^2 + S^2)))
 if(lower.tail==FALSE) result <- 1- result
 if(log.p==TRUE)       result <- log(result)
 return(result)
 }

