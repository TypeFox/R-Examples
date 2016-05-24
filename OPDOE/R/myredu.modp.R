# internal function from package crossdes
"myredu.modp" <-
function(a,b,p){

 la <- length(a)
 lb <- length(b)
 b1 <- b/b[1]                # divide by coefficient of highest power
  
 dummy <- numeric(la)        # partial result to be subtracted 
 
 for (i in 1:(la-lb+1)){
   dummy <- rep(0,la)
   dummy[i:(i+lb-1)] <- a[i]*b1
   a <- (a - dummy)%%p
 } 
 a[(la-lb+2):la]  
}
