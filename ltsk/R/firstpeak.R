firstpeak <-
function(dist,gamma)
{
 ## dist: distance intervals
 ## gamma : estimated semi-variogram ( with missing values)
 ## value: first peak of gamma
 ii <- which(is.na(gamma))
 jj <- which(!is.na(gamma))
 maxg <- max(gamma[jj])
 if( length(ii) > 0){
   if(length(jj)>1){
     f <- approxfun(dist, gamma,rule=2)
     gamma[ii] <- f(dist[ii])      
   }
   else{
     gamma[ii] <- maxg
   }
 }
 kk <- which( gamma > (maxg *.8))
 if(length(kk)>0){
   r <- min(kk)
 }
 else{
   r <- 1
 }
 r
}
