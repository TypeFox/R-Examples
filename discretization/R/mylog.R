mylog <-
function(x){
   x[which(x<=1.0e-10)] <- 1
   return (log(x))
}
