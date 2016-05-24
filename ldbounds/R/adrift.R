"adrift" <-
function(t,za,zb,pow){
   dr <- (zb[length(t)]+qnorm(pow))/sqrt(t[length(t)])
   drft <- bisect(t,za,zb,pow,dr)
   return(drft)
 }

