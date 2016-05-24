##
##  PURPOSE:   Square root of a general squared matrix
##
##  AUTHOR:    Arnost Komarek (LaTeX: Arno\v{s}t Kom\'arek)
##             arnost.komarek[AT]mff.cuni.cz
##
##  CREATED:   22/01/2008
##
##  FUNCTIONS:  MatSqrt
##
## ======================================================================

## *************************************************************
## MatSqrt
## *************************************************************
MatSqrt <- function(A)
{
  thispackage <- "mixAK"
  
  if (is.null(dim(A))){
    if (any(A < 0)) return(complex(real=sqrt(abs(A))*(A >= 0), imaginary=sqrt(abs(A))*(A < 0)))
    else            return(sqrt(A))
  }   
  
  p <- nrow(A)
  if (ncol(A) != p) stop("A must be a squared matrix")

  ldwork <- p*p

  RES <- .C("sqrtGE", Asqrt.re=as.double(A),    Asqrt.im=double(p*p),     Vinv.re=double(p*p),  Vinv.im=double(p*p),  complexRES=integer(1),
                      sqrt.lambda.re=double(p), sqrt.lambda.im=double(p), V.re=double(p*p),     V.im=double(p*p),
                      dwork=double(ldwork),     jpvt=integer(p),
                      err=as.integer(0),        p=as.integer(p),
            PACKAGE=thispackage)
  if (RES$err) stop("Matrix square root not computed.")
  
  if (RES$complexRES){
    Asqrt <- matrix(complex(real=RES$Asqrt.re, imaginary=RES$Asqrt.im), nrow=p, ncol=p)
    #V     <- matrix(complex(real=RES$V.re, imaginary=RES$V.im), nrow=p, ncol=p)    
    #Vinv  <- matrix(complex(real=RES$Vinv.re, imaginary=RES$Vinv.im), nrow=p, ncol=p)
    #slambda <- complex(real=RES$sqrt.lambda.re, imaginary=RES$sqrt.lambda.im)
  }else{
    Asqrt <- matrix(RES$Asqrt.re, nrow=p, ncol=p)
    #V     <- matrix(RES$V.re, nrow=p, ncol=p)        
    #Vinv  <- matrix(RES$Vinv.re, nrow=p, ncol=p)
    #slambda <- RES$sqrt.lambda.re
  }   
  
  rownames(Asqrt) <- rownames(A)
  colnames(Asqrt) <- colnames(A)

  return(Asqrt)                      
}  
