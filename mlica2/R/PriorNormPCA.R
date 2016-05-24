PriorNormPCA <-
function(X){ # START FUNCTION
 ndim <- ncol(X);
 ntp  <- nrow(X);
 # Center column means to zero
 for( s in 1:ndim ){
  X[,s] <- X[,s] - mean( X[,s] );
 }
 print("Performing SVD");
 # SVD
 svd.o <- svd(X,LINPACK=FALSE);
 Dx <- diag(svd.o$d*svd.o$d)/ntp;
 Ex <- svd.o$v ;
 barplot(Dx,main="Singular values");
 return(list(X=X,Dx=Dx,Ex=Ex));
}
