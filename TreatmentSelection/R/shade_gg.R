shade_gg <-
function(p, bounds, fixed, type, lty = 1, bands, width=5){
  
   low <- high <- y <- NULL #appease check
  
   if(any(dim(bounds)==0)){
     
     return(p )
   }
   
   
   yend <- xend <- NULL #to appease check
   
   numbounds <- ncol(bounds)
   if(is.null(numbounds)) {
      numbounds <- 1
      bounds <- matrix(bounds, ncol = 1)
   } 

   mydatERRORBAR <- data.frame(t(bounds), fixed)
   names(mydatERRORBAR) = c("low", "high", "fixed")

   mydat <- data.frame( c(bounds[1,], rev(bounds[2,])), c(fixed, rev(fixed)))
   names(mydat) = c("bounds", "fixed")
   if(type == "v"){
     #vertical ci
     
   if(bands == FALSE){
     p <- p + geom_errorbar(data = mydatERRORBAR, aes(ymin = low, ymax = high, x = fixed), linetype = lty, color = "gray50", size = 1, width = width)
   }else{
     #p <- p + geom_ribbon(data = mydat, alpha =.2, aes(x = fixed, ymin = low, ymax = high))
     p <- p + geom_polygon(data = mydat, alpha = .2, aes(y = bounds, x = fixed))
   }
   }else if(type =="h"){
   if(bands == FALSE){

      #p <- p + geom_errorbar(data = mydatERRORBAR, aes(ymin = low, ymax = high, x = fixed), linetype = lty, color = "gray50", size = 1, width = width) 
      p <- p + geom_segment(data = mydatERRORBAR, aes(y = fixed, yend = fixed, x = low, xend = high), linetype = lty, color = "gray50", size = 1)#, width = width)
      mydat$y = mydat$fixed-width/2; mydat$yend = mydat$fixed + width/2
      p <- p + geom_segment(data = mydat, aes(y = y, yend=yend, x = bounds, xend = bounds), linetype = lty, color = "gray50", size = 1)

      }else{

      #p <- p + geom_ribbon(data = mydat, alpha =.2, aes(x = fixed, ymin = low, ymax = high)) #+ scale_x_reverse()
      p <- p + geom_polygon(data = mydat, alpha = .2, aes(x = bounds, y = fixed))
   }
   }
   return(p)
   
}
