

##################
#### Gproduct ####
##################

 ## Description:
 ##     internal function in dbpls.
 ##         - calculate F for each iteration
 ##
 
 
Gproduct <- function(f,Dw,G0)
{
   f2 <- as.numeric(t(f)%*%Dw%*%f)
   g0 <- (G0%*%Dw%*%f)/f2
   g00 <- as.numeric((t(g0)%*%Dw%*%f)/f2)
   h0 <- g0-1/2*(g00*f)
   G1 <- G0-f%*%t(h0)-h0%*%t(f)
   return(G1)
}
 