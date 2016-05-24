`p4pl` <-
function(theta=0,a=1,b=0,c=0,d=1,lower.tail=TRUE,log.p=FALSE) {
 if(lower.tail==TRUE)  result <-  c + (d-c)/(1+exp(-a*(theta-b)))
 if(lower.tail==FALSE) result <- 1-(c + (d-c)/(1+exp(-a*(theta-b))))
 if(log.p==TRUE)       result <- log(result)
 return(result)
 }

