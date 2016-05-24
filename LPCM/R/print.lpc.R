print.lpc <-
  function(x, digits=max(3,getOption('digits')-3), ...){
      
      sx <- as.character(substitute(x))
      if (sum(nchar(sx))>200){ sx<-"name.of.this.object"}

      cat("\n")
      cat("Type plot(", sx, ") to see a graphical display of the fitted object. \n\n")
      cat("Type names(", sx, ") to see an overview of items available. \n\n")
              
      if(class(x)=="lpc"){
         if (x$scaled){
             cat("The data have been scaled by dividing through \n")
             cat(x$Misc$scaled.by)
             cat("\n")
         } else {
           cat("The data have not been scaled. \n")
         }
       }
     else if(class(x)=="lpc.spline"){
        cat("A cubic spline with ", dim(x$knots.coords[[1]])[2], " knots and total arc length ", diff(range(x$knots.pi[[1]])), " has been laid through the  local centers of mass representing the local principal curve. \n") 
     }
      
 }    



print.lpc.spline <-print.lpc
