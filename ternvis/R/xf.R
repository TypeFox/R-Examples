xf <-
function(p=cbind(1,1,1)/3,M=tsetup()$M32) {
   n <- nrow(p)
   xout <- NA * p[,1] # default for fault finding
   for (i in 1:n){    
    if(!is.na(p[i,2]))
     if(!is.na(p[i,3])) 
                  xout[i] <-  (M %*% p[i,])[1]
   }	    
   cbind(xout)
}
