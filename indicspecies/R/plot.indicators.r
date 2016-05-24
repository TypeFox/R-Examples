plot.indicators<-function(x, type="sqrtIV", maxline=TRUE, ...) {
   A = x$A
   B = x$B
   sqrtIV=x$sqrtIV
   order = rowSums(x$C)
   if(is.data.frame(A)) {
	   if(type=="IV") val = sqrtIV[,1]^2
   	   else if(type=="sqrtIV") val = sqrtIV[,1]
   	   else if(type=="A") val = A[,1]	
       else if(type=="B") val = B[,1]	
       else if(type=="LA") val = A[,2]	
       else if(type=="UA") val = A[,3]	
       else if(type=="LB") val = B[,2]	
       else if(type=="UB") val = B[,3]	
       else if(type=="LsqrtIV") val = sqrtIV[,2]	
       else if(type=="UsqrtIV") val = sqrtIV[,3]	
   } else {
	   if(type=="IV") val = sqrtIV^2
   	   else if(type=="sqrtIV") val = sqrtIV
   	   else if(type=="A") val = A	
       else if(type=="B") val = B	   	
   }
   plot(order,val, type="n", axes=FALSE, xlab="Order", ylab=type,...)	
   points(order,val, pch=1, cex=0.5)	
   axis(1, at = order, labels=order)
   axis(2)
   if(maxline) lines(1:max(order),tapply(val,order,max), col="gray")
}
