###############################################################################
hfun <- function(v, k) {
   
   result <- array(NA, dim=dim(v))
   smallers <- abs(v) < k
   
   result[smallers]  <- (v[smallers]^2)/2 * (1 - v[smallers]^2 / k^2 + v[smallers]^4 / (3*k^4) ) 
   result[!smallers] <- k^2 / 6
   
   return (result)
}
 
###############################################################################
phifun <- function (v, k) {

   result <- array(NA, dim=dim(v))
   smallers <- abs(v) < k
   
   result[smallers]  <- v[smallers] * ( 1-( v[smallers]^2 / k^2) )^2
   result[!smallers] <- 0
   
   return (result)
}
###############################################################################
dphifun <- function(v,k) {

   result <- array(NA, dim=dim(v))
   smallers <- abs(v) < k
   
   result[smallers]  <-  (sqrt(1-(v[smallers]^2/k^2))) - ( (v[smallers]^2/k^2) * (1-(v[smallers]^2/k^2)))
   result[!smallers] <-  0
   
   return (result)
}   
###############################################################################
#                         DELTA ERROR TAO
###############################################################################
deltaE.TAO <- function (arguments) {
   prediction <- arguments[[1]]
   target     <- arguments[[2]] 
   Stao       <- arguments[[3]]$deltaE$Stao  # the third argument is the net.
   
   residual   <- prediction - target
   scaled.residual <- residual / Stao
   c1  <- 1.56
   c2  <- 6.08
   bf  <- c1^2 / 12 

   h1  <- hfun(scaled.residual, c1)
   h2  <- hfun(scaled.residual, c2)
   phi1 <- phifun(scaled.residual,c1)
   phi2 <- phifun(scaled.residual,c2)
   
   if (sum(phi1 * residual) == 0.0) {
     dS2e <- 0.0
   } else {
     dS2e <- Stao * (sum(phi1) / (sum(phi1*residual)))
   }
   
   result <-mean(2*Stao*dS2e*h2 + phi2*(Stao - dS2e * residual))
   return(result)
   
}
###############################################################################
#                         DELTA ERROR LMS 
###############################################################################
deltaE.LMS <- function(arguments) {
   prediction <- arguments[[1]]                      # arg1 is the prediction
   target     <- arguments[[2]]                      # arg2 is the target
   residual   <- prediction - target
   return(residual)
}
###############################################################################
#                         DELTA ERROR LMLS
###############################################################################
deltaE.LMLS <- function(arguments) {
   prediction <- arguments[[1]]                      # arg1 is the prediction
   target     <- arguments[[2]]                      # arg2 is the target
   residual   <- prediction - target
   result     <- residual / (1 + residual^2 / 2) 
   return(result)
}
###############################################################################
#                         ERROR LMS 
###############################################################################
error.LMS <- function(arguments) {
   prediction <- arguments[[1]]                     # arg1 is the prediction
   target     <- arguments[[2]]                     # arg2 is the target
   residual   <- prediction - target
   result     <- mean((prediction - target)^2)
   return(result)
}
###############################################################################
#                         ERROR LMLS
###############################################################################
error.LMLS <- function(arguments) {
   prediction <- arguments[[1]]                     # arg1 is the prediction
   target     <- arguments[[2]]                     # arg2 is the target
   residual   <- prediction - target 
   result     <- mean(log(1 + residual^2 / 2))
   return(result)
}
###############################################################################
#                          ERROR TAO
###############################################################################
error.TAO <- function(arguments) {
   prediction <- arguments[[1]]                     # arg1 is the prediction
   target     <- arguments[[2]]                     # arg2 is the target
   Stao       <- arguments[[3]]$deltaE$Stao # arg3 is net
   residual   <- prediction - target 
   
   n.residual <- nrow(residual)
   perf <- NA
   
   scaled.residual <- residual / Stao
   c1  <- 1.56
   c2  <- 6.08
   bf  <- c1^2 / 12 

   h1  <- hfun(scaled.residual, c1)
   h2  <- hfun(scaled.residual, c2)
   
   
   new.Stao <- Stao*sqrt(sum(h1)/(n.residual * bf))  # n.residuals o n.residuals*n.output.MLPneurons ??
   tao.error.squared <- new.Stao^2 * mean(h2)
   return(list(perf=tao.error.squared, Stao=new.Stao))
}
###############################################################################

    
            
	    


