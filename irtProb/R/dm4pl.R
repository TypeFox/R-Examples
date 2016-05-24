`dm4pl` <-
function(theta=0,S=rep(0,length(theta)),C=rep(0,length(theta)),D=rep(0,length(theta)),b=0,s=1/1.702,c=0,d=1,log.p=FALSE) {
 derpm4pl  <- deriv3(dens ~ (C+c) + ((d-D) - (C+c))/(1+exp(-(theta-b)/sqrt(s^2 + S^2))) , c("theta","S","C","D","s","b","c","d"), func = TRUE, hessian = FALSE) 
 result    <- attr(derpm4pl(theta=theta,S=S,C=C,D=D,s=s,b=b,c=c,d=d), "gradient")[,1]
 if (log.p==TRUE) result <- log(result)
 return(result)
 }

