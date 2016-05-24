"parcau" <-
function(lmom) {
   para <- vector(mode="numeric", length=2)
   names(para) <- c("xi","alpha")
   if(length(lmom$source) == 1 && lmom$source != "TLmoms" ) {
     warning("TL-moments with trim=1 are required--can not complete parameter estimation")
     return()    
   }
   if(lmom$trim != 1) {
     warning("Attribute of TL-moments is not trim=1--can not complete parameter estimation")
     return()      
   }
   para[1] <- lmom$lambdas[1] 
   para[2] <- lmom$lambdas[2]/0.698 
   return(list(type = 'cau', para=para, source="parcau"))
}
