print.summary.tmatrix <-
function(x,...){
  
  
   cat(paste("transition matrix",x$name.mat,"\n\n"))
   cat(paste("lambda:", round(x$lambda, 3),"\n\n"))
   cat("stable stage distribution: \n")
   names(x$stable.stage.distribution)<- x$m.names
   print(noquote(format(round(x$stable.stage.distribution,3)))); cat("\n\n")
   
   
   cat("reproductive value:\n")
   names(x$reproductive.value)<- x$m.names
   print(noquote(format(round(x$reproductive.value,3)))); cat("\n\n")
    
   cat("sensitivities:\n")
   dimnames(x$sensitivity) <-list(x$m.names, x$m.names)
   print(round(x$sensitivity,3)); cat("\n\n")
 
  cat("elasticities:\n")
  class(x$elasticity) <-"matrix"
  dimnames(x$elasticity) <-list(x$m.names, x$m.names)
  print(round(x$elasticity,3)); cat("\n\n")
    
 }

