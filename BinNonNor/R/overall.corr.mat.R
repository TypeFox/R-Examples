overall.corr.mat <-
function(n.BB, n.NN, prop.vec=NULL, corr.vec = NULL, corr.mat=NULL, coef.mat=NULL){

   validation.bin(n.BB, prop.vec)

   if(is.null(corr.mat) && !is.null(corr.vec)) {
   d=ceiling(uniroot(function(d) d^2-d-2*length(corr.vec), interval=c(0,1000))$root)
   corr.mat=diag(1,d)
   corr.mat[lower.tri(corr.mat)]=corr.vec
   corr.mat=corr.mat+t(corr.mat)-diag(1,d)
   }

   correlation.bound.check(n.BB, n.NN, prop.vec, corr.vec=NULL, corr.mat, coef.mat)

   if(!is.null(prop.vec) && is.null(coef.mat) ) {
   final.corr.mat<-diag(1,n.BB)
   final.corr.mat[1:n.BB,1:n.BB]=Tetra.Corr.BB(n.BB, prop.vec, corr.vec = NULL, corr.mat)
   } else
   if( is.null(prop.vec) && !is.null(coef.mat) ) {
   final.corr.mat<-diag(1,n.NN)
   final.corr.mat[1:n.NN,1:n.NN]=Int.Corr.NN(n.NN, corr.vec = NULL, corr.mat, coef.mat) 
   } else
   if(!is.null(prop.vec) && !is.null(coef.mat) ) {
   final.corr.mat<-diag(1,(n.BB+n.NN))
   final.corr.mat[1:n.BB,1:n.BB]=Tetra.Corr.BB(n.BB, prop.vec, corr.vec = NULL, corr.mat)
   final.corr.mat[(n.BB+1):(n.BB+n.NN),(n.BB+1):(n.BB+n.NN)]=Int.Corr.NN(n.NN, corr.vec = NULL, corr.mat, coef.mat)
   final.corr.mat[1:n.BB,(n.BB+1):(n.BB+n.NN)]=Biserial.Corr.BN(n.BB, n.NN, prop.vec, corr.vec = NULL, corr.mat, coef.mat)
   final.corr.mat[(n.BB+1):(n.BB+n.NN),1:n.BB]=t(Biserial.Corr.BN(n.BB, n.NN, prop.vec, corr.vec = NULL, corr.mat, coef.mat))
   }

   if(is.positive.definite(final.corr.mat)==FALSE) {
     warning("Intermediate correlation matrix is not positive definite. Nearest positive definite matrix is used!")
     final.cor.mat = as.matrix(nearPD(final.corr.mat, corr = TRUE, keepDiag =TRUE)$mat)
   }

return(final.corr.mat)
}
