lpc.project <- function(object, newdata, ...){
  if (class(object)=="lpc"){
       lpcobject  <- object
       lpcsl  <- lpc.splinefun(lpcobject)
   } else {
       lpcobject  <- object$lpcobject
       lpcsl  <- object$splinefun
   }

   if (lpcobject$scaled){
      data <-sweep(as.matrix(newdata),2, lpcobject$Misc$scaled.by, "/")
   }   else {
      data <-newdata
    }
  
  
   result <- lpc.project.spline(lpcsl, data,...)[c(2,3,4,5,1)]
  result
}    
    
    
    

    
