`disattenuated.cor` <-
function(r.xy, r.xx, new.r.xx=1){

   new.r <- function(rxy, rxx, nrxx) sqrt(prod(nrxx))*rxy/sqrt(prod(rxx))
   
   if(! is.matrix(r.xy)){
     new.r.xy <- new.r(r.xy, r.xx, new.r.xx)  
   } else{
     new.r.xy <- r.xy
     for(i in 2:ncol(r.xy)){
      for(j in 1:(i-1)){
       if(length(new.r.xx)==1) nr <- new.r.xx else nr <- new.r.xx[c(i,j)]
       new.r.xy[i,j] <- new.r(r.xy[i,j], r.xx[c(i,j)], nr) 
      }
     }
     diag(new.r.xy) <- r.xx
     }
   t(new.r.xy)
}

