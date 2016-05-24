correlation.bound.check <-
function(n.P, n.B, n.O, lambda.vec=NULL, prop.vec=NULL, prop.list=NULL,corr.vec = NULL, corr.mat = NULL) {
 
   if(is.null(corr.mat) && !is.null(corr.vec)) {
   d = n.P+n.B+n.O
   corr.mat=diag(1,d)
   corr.mat[lower.tri(corr.mat)]=corr.vec
   corr.mat=corr.mat+t(corr.mat)-diag(1,d)
   }

   validation.corr(n.P, n.B, n.O, corr.vec=NULL, corr.mat)
   limitscor.mat=correlation.limits(n.P,n.B,n.O, lambda.vec, prop.vec, prop.list)
 
   if((n.P+n.B+n.O)>1) {
   errorCount= 0
   for (i in 1:(n.P+n.B+n.O-1))   {
   for (j in (i+1):(n.P+n.B+n.O)) {
   if (i != j) {
   if (corr.mat[i, j] > limitscor.mat[i, j] | corr.mat[j, i] <  limitscor.mat[j, i]) {
         cat("\n corr.mat[", i, ",", j, "] must be between ", round(limitscor.mat[j,i], 7), " and ", round(limitscor.mat[i,j], 7), "\n")
         errorCount = errorCount + 1
         cat("\n")
   }
   }
   }
   }
   
   if (errorCount > 0) {
   stop("Range violation occurred in the target correlation matrix! \n")
   }
   }
return(TRUE)
}
