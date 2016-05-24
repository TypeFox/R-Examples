E.mat <- function(X,min.MAF=NULL,max.missing=NULL,impute.method="mean",tol=0.02,
                  n.core=1,shrink=FALSE,return.imputed=FALSE, type="A#A", ploidy=2){
  
  if(type == "A#A"){
    X4 <- A.mat(X, min.MAF=min.MAF,max.missing=max.missing,impute.method=impute.method,tol=tol,
                n.core=n.core,shrink=shrink,return.imputed=return.imputed, ploidy=ploidy)
    X5 <- matrixcalc::hadamard.prod(X4,X4)
  }
  if(type == "A#D"){
    X4 <- A.mat(X, min.MAF=min.MAF,max.missing=max.missing,impute.method=impute.method,tol=tol,
                n.core=n.core,shrink=shrink,return.imputed=return.imputed, ploidy=ploidy)
    X4D <- D.mat(X, min.MAF=min.MAF,max.missing=max.missing,impute.method=impute.method,tol=tol,
                       n.core=n.core,shrink=shrink,return.imputed=return.imputed,ploidy=ploidy)
    X5 <- matrixcalc::hadamard.prod(X4,X4D)
  }
  if(type == "D#D"){
    X4D <- D.mat(X, min.MAF=min.MAF,max.missing=max.missing,impute.method=impute.method,tol=tol,
                 n.core=n.core,shrink=shrink,return.imputed=return.imputed, ploidy=ploidy)
    X5 <- matrixcalc::hadamard.prod(X4D,X4D)
  }
  return(X5)
}