orank.default <-
function(veg,use,rlimit=5,y=1,x.axis=NULL,y.axis=NULL,...) {
#  revised with new default, 25. 6. 2014
   o.orank<- orank1(veg,use,rlimit,y,x.axis,y.axis)
   class(o.orank)<- "orank"
   o.orank$call<- match.call()
   cat("Call:\n") 
   print(o.orank$call)
   cat("\n")
   cat("RANK","\n")
   cat("Variable              rank no.    var.      var.%   cum. var.%","\n")
  for(i in 1:o.orank$n.ranks){
     cat(sprintf("%-20s %4i %12.2f %10.3f %12.3f",o.orank$var.names[i],i,o.orank$var.explained[i],o.orank$var.percent[i],o.orank$cum.var[i]))
     cat("\n")
  }
   o.orank
 }
