srank.default <-
function(veg,groups,method,y,...) {
     o.srank <- srank2(veg,groups,method,y)
     class(o.srank)<- "srank"
     o.srank$call<- match.call()
     cat("Call:\n") 
     print(o.srank$call)
     o.srank
 }
