summary.orank <-
function(object,...) {
   o.orank<- object
   cat("RANK","\n")
   cat("Variable              rank no.    var.      var.%   cum. var.%","\n")
  for(i in 1:o.orank$n.ranks){
     cat(sprintf("%-20s %4i %12.2f %10.3f %12.3f",o.orank$var.names[i],i,o.orank$var.explained[i],o.orank$var.percent[i],o.orank$cum.var[i]))
     cat("\n")
  }
}
