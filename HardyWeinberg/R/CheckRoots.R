`CheckRoots` <- function(roots,limits=c(0,1),mini=TRUE,verbose=FALSE) {
if(verbose) {  
   cat("checking roots\n")
   print(roots)
   }
   roots <- roots[round(roots$rAB,digits=5)==round(roots$check,digits=5),]
   roots <- roots[roots$"Im<1e-10"==TRUE,]
   roots <- roots[(roots$Re < limits[2] & roots$Re > limits[1]),]
   roots <- roots$Re
   if(length(roots)>1)
     if(verbose) cat("multiple roots\n")
   if(length(roots)!= 0) {
     if(mini)
        y <- min(roots)
     else y <- max(roots)
   }
   else y<-NULL
   if(verbose) print(y)
   return(y)
}

