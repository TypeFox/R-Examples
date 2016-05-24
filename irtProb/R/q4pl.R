`q4pl` <-
function(p=0.05,a=1,b=0,c=0,d=1,lower.tail=TRUE,log.p=FALSE) {
 p      <- ifelse(log.p==FALSE, p, exp(p))
 p      <- ifelse(lower.tail==TRUE, p, 1-p)
 result <- log((d-c)/(p-c) - 1)/(-a) + b
 return(result)
 }

