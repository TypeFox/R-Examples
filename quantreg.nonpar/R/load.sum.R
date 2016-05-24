load.sum <-
function(vec){
  
 if (is.factor(vec)){
   if(is.ordered(vec)){
     return (median(vec))
   } else {
     return (names(sort(-table(vec)))[1])
   }
 } else {
   return (mean(vec))
 }
}
