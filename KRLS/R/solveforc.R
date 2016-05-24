solveforc <-
function(y=NULL,Eigenobject=NULL,lambda=NULL,eigtrunc=NULL){
	
	nn <- nrow(y)
	Ginv <- matrix(NA,nn,nn)
	
	if (is.null(eigtrunc)){	Ginv <- tcrossprod(multdiag(X=Eigenobject$vectors,d=1/(Eigenobject$values+lambda)),Eigenobject$vectors)
} else{
	
	#eigentruncation: keep only eigenvectors at least 'eigtrunc' times as large as the largest
	lastkeeper=max(which(Eigenobject$values>=eigtrunc*Eigenobject$values[1]))
  #We need at least one!
  lastkeeper=max(1,lastkeeper)
  
	Ginv <- tcrossprod(multdiag(X=Eigenobject$vectors[,1:lastkeeper],d=1/(Eigenobject$values[1:lastkeeper]+lambda)),Eigenobject$vectors[,1:lastkeeper])
}
		
	coeffs <- matrix(NA,nrow=nn,ncol=1)
  coeffs <- tcrossprod(Ginv,t(y))
	Le      <- crossprod(coeffs/diag(Ginv))
	#Le      <- as.vector(crossprod(coeffs/diag(Ginv)))
  return(list(coeffs=coeffs,
               Le=Le))
  }


