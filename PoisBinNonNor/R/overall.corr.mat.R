overall.corr.mat <-
function(n.P, n.B, n.C, lambda.vec=NULL, prop.vec=NULL, coef.mat=NULL, corr.vec = NULL, corr.mat=NULL){

   validation.bin(n.B, prop.vec)
 
   if(is.null(corr.mat) && !is.null(corr.vec)) {
   d=ceiling(uniroot(function(d) d^2-d-2*length(corr.vec), interval=c(0,1000))$root)
   corr.mat=diag(1,d)
   corr.mat[lower.tri(corr.mat)]=corr.vec
   corr.mat=corr.mat+t(corr.mat)-diag(1,d)
   }

   correlation.bound.check(n.P, n.B, n.C, lambda.vec, prop.vec, coef.mat, corr.vec=NULL, corr.mat=corr.mat)

   if(!is.null(lambda.vec) && is.null(prop.vec) && is.null(coef.mat) ) {
   final.corr.mat<-diag(1,n.P)
   final.corr.mat[1:n.P,1:n.P]=intermediate.corr.PP(n.P, n.B, n.C, lambda.vec, corr.vec=NULL, corr.mat)
   } else
   if(is.null(lambda.vec) && !is.null(prop.vec) && is.null(coef.mat) ) {
   final.corr.mat<-diag(1,n.B)
   final.corr.mat[1:n.B,1:n.B]=intermediate.corr.BB(n.P,n.B,n.C, prop.vec, corr.vec = NULL, corr.mat)
   } else
   if(is.null(lambda.vec) && is.null(prop.vec) && !is.null(coef.mat) ) {
   final.corr.mat<-diag(1,n.C)
   final.corr.mat[1:n.C,1:n.C]=intermediate.corr.CC(n.P, n.B, n.C, coef.mat, corr.vec = NULL, corr.mat)
   } else
   if(!is.null(lambda.vec) && !is.null(prop.vec) && is.null(coef.mat) ){
   final.corr.mat<-diag(1,(n.P+n.B))
   final.corr.mat[1:n.P,1:n.P]=intermediate.corr.PP(n.P, n.B, n.C, lambda.vec, corr.vec=NULL, corr.mat)
   final.corr.mat[(n.P+1):(n.P+n.B),(n.P+1):(n.P+n.B)]=intermediate.corr.BB(n.P,n.B,n.C, prop.vec, corr.vec = NULL, corr.mat)
   final.corr.mat[1:n.P,(n.P+1):(n.P+n.B)]=intermediate.corr.PB(n.P, n.B, n.C, lambda.vec, prop.vec, coef.mat, corr.vec= NULL, corr.mat)
   final.corr.mat[(n.P+1):(n.P+n.B),1:n.P]=t(final.corr.mat[1:n.P,(n.P+1):(n.P+n.B)])
   } else
   if(!is.null(lambda.vec) && is.null(prop.vec) && !is.null(coef.mat) ) {
   final.corr.mat<-diag(1,(n.P+n.C))
   final.corr.mat[1:n.P,1:n.P]=intermediate.corr.PP(n.P, n.B, n.C, lambda.vec, corr.vec=NULL, corr.mat)
   final.corr.mat[(n.P+1):(n.P+n.C),(n.P+1):(n.P+n.C)]=intermediate.corr.CC(n.P, n.B, n.C, coef.mat, corr.vec = NULL, corr.mat)
   final.corr.mat[1:n.P,(n.P+1):(n.P+n.C)]=intermediate.corr.PC(n.P, n.B, n.C, lambda.vec, prop.vec, coef.mat, corr.vec = NULL, corr.mat)
   final.corr.mat[(n.P+1):(n.P+n.C),1:n.P]=t(final.corr.mat[1:n.P,(n.P+1):(n.P+n.C)])
   } else
   if(is.null(lambda.vec) && !is.null(prop.vec) && !is.null(coef.mat) ) {
   final.corr.mat<-diag(1,(n.B+n.C))
   final.corr.mat[1:n.B,1:n.B]=intermediate.corr.BB(n.P,n.B,n.C, prop.vec, corr.vec = NULL, corr.mat)
   final.corr.mat[(n.B+1):(n.B+n.C),(n.B+1):(n.B+n.C)]=intermediate.corr.CC(n.P, n.B, n.C, coef.mat, corr.vec = NULL, corr.mat)
   final.corr.mat[1:n.B,(n.B+1):(n.B+n.C)]=intermediate.corr.BC(n.P, n.B, n.C, lambda.vec, prop.vec, coef.mat, corr.vec = NULL, corr.mat)
   final.corr.mat[(n.B+1):(n.B+n.C),1:n.B]=t(final.corr.mat[1:n.B,(n.B+1):(n.B+n.C)])
   } else
   if(!is.null(lambda.vec) && !is.null(prop.vec) && !is.null(coef.mat) ) {
   final.corr.mat<-diag(1,(n.P+n.B+n.C))
   final.corr.mat[1:n.P,1:n.P]=intermediate.corr.PP(n.P, n.B, n.C, lambda.vec, corr.vec=NULL, corr.mat)
   final.corr.mat[(n.P+1):(n.P+n.B),(n.P+1):(n.P+n.B)]=intermediate.corr.BB(n.P,n.B,n.C, prop.vec, corr.vec = NULL, corr.mat)
   final.corr.mat[(n.P+n.B+1):(n.P+n.B+n.C),(n.P+n.B+1):(n.P+n.B+n.C)]=intermediate.corr.CC(n.P, n.B, n.C, coef.mat, corr.vec = NULL, corr.mat)
   final.corr.mat[1:n.P,(n.P+1):(n.P+n.B)]=intermediate.corr.PB(n.P, n.B, n.C, lambda.vec, prop.vec, coef.mat, corr.vec= NULL, corr.mat)
   final.corr.mat[(n.P+1):(n.P+n.B),1:n.P]=t(final.corr.mat[1:n.P,(n.P+1):(n.P+n.B)])
   final.corr.mat[1:n.P,(n.P+n.B+1):(n.P+n.B+n.C)]=intermediate.corr.PC(n.P, n.B, n.C, lambda.vec, prop.vec, coef.mat, corr.vec = NULL, corr.mat)
   final.corr.mat[(n.P+n.B+1):(n.P+n.B+n.C),1:n.P]=t(final.corr.mat[1:n.P,(n.P+n.B+1):(n.P+n.B+n.C)])
   final.corr.mat[(n.P+1):(n.P+n.B),(n.P+n.B+1):(n.P+n.B+n.C)]=intermediate.corr.BC(n.P, n.B, n.C, lambda.vec, prop.vec, coef.mat, corr.vec = NULL, corr.mat)
   final.corr.mat[(n.P+n.B+1):(n.P+n.B+n.C),(n.P+1):(n.P+n.B)]=t(final.corr.mat[(n.P+1):(n.P+n.B),(n.P+n.B+1):(n.P+n.B+n.C)])
   }
  
   if(is.positive.definite(final.corr.mat)==FALSE) {
     warning("Intermediate correlation matrix is not positive definite. Nearest positive definite matrix is used!")
     final.corr.mat = as.matrix(nearPD(final.corr.mat, corr = TRUE, keepDiag =TRUE)$mat)
   }

return(final.corr.mat)
}
