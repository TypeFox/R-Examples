`d4pl` <-
function(theta=0,a=1,b=0,c=0,d=1,log.p=FALSE) {
 derp4pl  <- deriv3(dens ~ c + (d-c)/(1+exp(-a*(theta-b))), c("theta","a","b","c","d"), func = TRUE, hessian = TRUE); res <- derp4pl(theta=theta,a=a,b=b,c=c,d=d)
 result   <- attr(res,"gradient")[,1]
 # result   <- p4pl(theta=theta,a=a,b=b,c=c,d=d)*(1-p4pl(theta=theta,a=a,b=b,c=c,d=d))
 if (log.p==TRUE) result <- log(result)
 return(result)
 }

