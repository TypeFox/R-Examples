correlation.bound.check <-
function(n.BB, n.NN, prop.vec=NULL, corr.vec = NULL, corr.mat = NULL, coef.mat=NULL) {
 
   validation.corr(n.BB, n.NN, corr.vec, corr.mat) 

   limitscor.mat=correlation.limits(n.BB,n.NN,prop.vec,coef.mat)

   if(is.null(corr.mat) && !is.null(corr.vec)) {
   d = n.BB + n.NN
   corr.mat=diag(1,d)
   corr.mat[lower.tri(corr.mat)]=corr.vec
   corr.mat=corr.mat+t(corr.mat)-diag(1,d)
   }

   if((n.BB+n.NN)>1) {
   errorCount= 0

   for (i in 1:(n.BB+n.NN-1))   {
   for (j in (i+1):(n.BB+n.NN)) {
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
