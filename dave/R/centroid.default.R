centroid.default <-
function(nveg,grel,y,...) {
   o.centroid<- rcentroid(nveg,grel,y,...)
   class(o.centroid)<- "centroid"
   o.centroid$call<- match.call()
   cat("Call:\n") 
   print(o.centroid$call)
   o.centroid
 }
