overall.corr.mat <-
function(n.P, n.B, n.O, lambda.vec=NULL, prop.vec=NULL, prop.list=NULL, corr.vec = NULL, corr.mat=NULL){

   validation.bin(n.B, prop.vec)
   validation.ord(n.O, prop.list)

   if(is.null(corr.mat) && !is.null(corr.vec)) {
   d=ceiling(uniroot(function(d) d^2-d-2*length(corr.vec), interval=c(0,1000))$root)
   corr.mat=diag(1,d)
   corr.mat[lower.tri(corr.mat)]=corr.vec
   corr.mat=corr.mat+t(corr.mat)-diag(1,d)
   }

   correlation.bound.check(n.P, n.B, n.O, lambda.vec, prop.vec, prop.list, corr.vec= NULL, corr.mat)

   if(!is.null(lambda.vec) && is.null(prop.vec) && is.null(prop.list) ) {
   final.corr.mat<-diag(1,n.P)
   final.corr.mat[1:n.P,1:n.P]=intermediate.corr.PP(n.P, lambda.vec, corr.vec=NULL, corr.mat)
   } else
   if(is.null(lambda.vec) && !is.null(prop.vec) && is.null(prop.list) ) {
   final.corr.mat<-diag(1,n.B)
   final.corr.mat[1:n.B,1:n.B]=intermediate.corr.BO(n.B, n.O=0, prop.vec, prop.list=NULL, corr.vec = NULL, corr.mat)
   } else
   if(is.null(lambda.vec) && is.null(prop.vec) && !is.null(prop.list) ) {
   final.corr.mat<-diag(1,n.O)
   final.corr.mat[1:n.O,1:n.O]=intermediate.corr.BO(n.B=0, n.O, prop.vec=NULL, prop.list, corr.vec = NULL, corr.mat)
   } else
   if(!is.null(lambda.vec) && !is.null(prop.vec) && is.null(prop.list) ) {
   final.corr.mat<-diag(1,(n.P+n.B))
   final.corr.mat[1:n.P,1:n.P]=intermediate.corr.PP(n.P, lambda.vec, corr.vec=NULL, corr.mat)
   final.corr.mat[(n.P+1):(n.P+n.B),(n.P+1):(n.P+n.B)]=intermediate.corr.BO(n.B, n.O=0, prop.vec, prop.list=NULL, corr.vec = NULL, corr.mat)
   final.corr.mat[1:n.P,(n.P+1):(n.P+n.B)]=intermediate.corr.P_BO(n.P, n.B, n.O=0, lambda.vec, prop.vec, prop.list=NULL, corr.vec = NULL, corr.mat)
   final.corr.mat[(n.P+1):(n.P+n.B),1:n.P]=t(intermediate.corr.P_BO(n.P, n.B, n.O=0, lambda.vec, prop.vec, prop.list=NULL, corr.vec = NULL, corr.mat))
   } else
   if(!is.null(lambda.vec) && is.null(prop.vec) && !is.null(prop.list) ) {
   final.corr.mat<-diag(1,(n.P+n.O))
   final.corr.mat[1:n.P,1:n.P]=intermediate.corr.PP(n.P, lambda.vec, corr.vec=NULL, corr.mat)
   final.corr.mat[(n.P+1):(n.P+n.O),(n.P+1):(n.P+n.O)]=intermediate.corr.BO(n.B=0, n.O, prop.vec=NULL, prop.list, corr.vec = NULL, corr.mat)
   final.corr.mat[1:n.P,(n.P+1):(n.P+n.O)]=intermediate.corr.P_BO(n.P, n.B=0, n.O, lambda.vec, prop.vec=NULL, prop.list, corr.vec = NULL, corr.mat)
   final.corr.mat[(n.P+1):(n.P+n.O),1:n.P]=t(intermediate.corr.P_BO(n.P, n.B=0, n.O, lambda.vec, prop.vec=NULL, prop.list, corr.vec = NULL, corr.mat))
   } else
   if(is.null(lambda.vec) && !is.null(prop.vec) && !is.null(prop.list) ) {
   final.corr.mat<-diag(1,(n.B+n.O))
   final.corr.mat[1:(n.B+n.O),1:(n.B+n.O)]=intermediate.corr.BO(n.B, n.O, prop.vec, prop.list, corr.vec = NULL, corr.mat)
   } else
   if(!is.null(lambda.vec) && !is.null(prop.vec) && !is.null(prop.list) ) {
   final.corr.mat<-diag(1,(n.P+n.B+n.O))
   final.corr.mat[1:n.P,1:n.P]=intermediate.corr.PP(n.P, lambda.vec, corr.vec=NULL, corr.mat)
   final.corr.mat[(n.P+1):(n.P+n.B+n.O),(n.P+1):(n.P+n.B+n.O)]=intermediate.corr.BO(n.B, n.O, prop.vec, prop.list, corr.vec = NULL, corr.mat)
   final.corr.mat[1:n.P,(n.P+1):(n.P+n.B+n.O)]=intermediate.corr.P_BO(n.P, n.B, n.O, lambda.vec, prop.vec, prop.list, corr.vec = NULL, corr.mat)
   final.corr.mat[(n.P+1):(n.P+n.B+n.O),1:n.P]=t(intermediate.corr.P_BO(n.P, n.B, n.O, lambda.vec, prop.vec, prop.list, corr.vec = NULL, corr.mat))
   }
  
   if(is.positive.definite(final.corr.mat)==FALSE) {
     warning("Intermediate correlation matrix is not positive definite. Nearest positive definite matrix is used!")
     final.corr.mat = as.matrix(nearPD(final.corr.mat, corr = TRUE, keepDiag =TRUE)$mat)
   }

return(final.corr.mat)
}
