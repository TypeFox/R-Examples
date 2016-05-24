shade <-
function(bounds, fixed, type, lty = 1, bands){

   numbounds <- ncol(bounds)
   if(is.null(numbounds)) {
      numbounds <- 1
      bounds <- matrix(bounds, ncol = 1)
   } 

   if(type == "v"){
   if(bands == FALSE){
     for( i in 1:numbounds) { 
       lines( rep(fixed[i], 2), bounds[1:2, i], lwd = 1.5, col = "gray50", lty = lty) 
       points( rep(fixed[i], 2), bounds[1:2, i], pch = "_", col = "gray50", cex = 1.1)
     }
   }else{
     
     polygon( c(fixed, rev(fixed)), c(bounds[1,], rev(bounds[2,])),density = 25, angle = 90, col="grey50", border = NA, lty = lty  )
     lines( fixed, bounds[1,], col="grey50") 
     lines( fixed, bounds[2,], col="grey50")
   
   }
   }else if(type =="h"){
   if(bands == FALSE){
      for(i in 1:numbounds) { 
         lines( bounds[1:2, i], rep(fixed[i], 2), lwd = 1.5,col = "gray50", lty = lty) 
         points( bounds[1:2, i], rep(fixed[i], 2), pch = "|", col = "gray50", cex = 1.1)
      }
   }else{
      polygon( c(bounds[1,], rev(bounds[2,])), c(fixed,rev(fixed)),density = 25, angle = 0, col="grey50", border = NA , lty = lty  )
      lines( bounds[1,], fixed,  col="grey") 
      lines( bounds[2,], fixed,  col="grey")
   
      
   }
   }
}
