#' Optimization & predictive modelling Toolsets
#' 
#' @description optR function for solving linear systems using numerical approaches. 
#' Current toolbox supports Gauss Elimination, LU decomposition, Conjugate Gradiant Decent and Gauss-Sideal methods for solving the system of form AX=b
#' For optimization using numerical methods cgm method performed faster in comparision with gaussseidel.
#' For decomposition LU is utilized for multiple responses to enhance the speed of computation.
#' @param x   : Input matrix
#' @param ... : S3 method
#' @return optR   : Return optR class
#' @author PKS Prakash
#' @export
#' @examples
#' # Solving equation Ax=b
#' A<-matrix(c(6,-4,1, -4,6,-4,1,-4,6), nrow=3,ncol=3, byrow = TRUE)
#' b<-matrix(c(-14,36, 6), nrow=3,ncol=1,byrow=TRUE)
#' Z<-optR(A, b, method="gauss") # Solve Linear model using Gauss Elimination
#' 
#' # Solve Linear model using LU decomposition (Supports Multi-response)
#' Z<-optR(A, b, method="LU") 
#' 
#' # Solve the matrix using Gauss Elimination (1, -1, 2)
#' A<-matrix(c(2,-2,6, -2,4,3,-1,8,4), nrow=3,ncol=3, byrow = TRUE)
#' b<-matrix(c(16,0, -1), nrow=3,ncol=1,byrow=TRUE)
#' Z<-optR(A, b, method="gauss") # Solve Linear model using Gauss Elimination
#' 
#' require(utils)
#' set.seed(129)
#' n <- 10 ; p <- 4
#' X <- matrix(rnorm(n * p), n, p) # no intercept!
#' y <- rnorm(n)
#' Z<-optR(X, y, method="cgm")
optR<-function(x, ...) UseMethod("optR")


#' Optimization & predictive modelling Toolsets 
#' 
#' @description optR package to perform the optimization using numerical methods
#' @param formula   : formula to build model
#' @param data      : data used to build model
#' @param method    : "gauss" for gaussian elimination and "LU" for LU factorization
#' @param iter  : Number of Iterations
#' @param tol   : Convergence tolerance 
#' @param ...   : S3 Class
#' @return U        : Decomposed matrix for Gauss-ELimination Ax=b is converted into Ux=c where U is upper triangular matrix for LU decomposition U contain the values for L & U decomposition LUx=b
#' @return c        : transformed b & for LU transformation c is y from equation Ux=y
#' @return estimates  : Return x values for linear system
#' @author PKS Prakash
#' @export
#' @examples
#' # Solving equation Ax=b
#' b<-matrix(c(-14,36, 6), nrow=3,ncol=1,byrow=TRUE)
#' A<-matrix(c(6,-4,1, -4,6,-4,1,-4,6), nrow=3,ncol=3, byrow = TRUE)
#' Z<-optR(b~A-1, method="gauss") # -1 to remove the constant vector
#' 
#' Z<-optR(b~A-1, method="LU") # -1 to remove the constant vector
#' 
#' require(utils)
#' set.seed(129)
#' n <- 10 ; p <- 4
#' X <- matrix(rnorm(n * p), n, p) # no intercept!
#' y <- rnorm(n)
#' data<-cbind(X, y)
#' colnames(data)<-c("var1", "var2", "var3", "var4", "y")
#' Z<-optR(y~var1+var2+var3+var4+var1*var2-1, data=data.frame(data), method="cgm")
optR.formula<-function(formula, data=list(), method=c("gauss, LU, gaussseidel", "cgm", "choleski"), iter=500, tol=1e-7, ...)
{
  # Parse the call
  cl <- match.call()
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "subset", "na.action"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame())
  
  # Extract data
  mt <- attr(mf, "terms")
  x <- model.matrix(mt, mf, contrasts)
  y<-model.response(mf, "numeric")
  
  # Default Method
  if(length(method)>1) method="cgm"
  
  # Fit Models
  nROWx=nrow(x)
  nCOLx=ncol(x)
  if(nROWx==nCOLx){
    optR<-optR.fit(x, y, method, iter, tol)  # Fit optimization method
  } else {
    y<-t(x)%*%y
    x<-t(x)%*%x
    optR<-optR.fit(x, y, method, iter, tol) # Fit optimization method
  }
  
  
  optR$formula<-formula
  optR$na.action <- attr(mf, "na.action")
  optR$xlevels <- .getXlevels(mt, mf)
  optR$terms <- mt
  optR$call<-cl
  optR$method<-method
  class(optR)<-"optR"
  optR
  # Call the default function
}


#' Optimization & predictive modelling Toolsets
#' 
#' optR is the default function for optimization
#' @param x     : Input data frame
#' @param y     : Response is data frame
#' @param method  : "gauss" for gaussian elimination and "LU" for LU factorization
#' @param iter  : Number of Iterations
#' @param tol   : Convergence tolerance
#' @param ...   : S3 Class
#' @return U    : Decomposed matrix for Gauss-ELimination Ax=b is converted into Ux=c where U is upper triangular matrix for LU decomposition U contain the values for L & U decomposition LUx=b   
#' @return c    : transformed b & for LU transformation c is y from equation Ux=y
#' @return estimates  : Return x values for linear system
#' @return seq        : sequence of A matrix re-ordered
#' @export
#' @examples
#' # Solving equation Ax=b
#' A<-matrix(c(6,-4,1, -4,6,-4,1,-4,6), nrow=3,ncol=3, byrow = TRUE)
#' b<-matrix(c(-14,36, 6), nrow=3,ncol=1,byrow=TRUE)
#' Z<-optR(A, b, method="gauss") 
#' 
#' # Solve Linear model using LU decomposition (Supports Multi-response)
#' Z<-optR(A, b, method="LU")
#' 
#' # Solving the function using numerical method
#' Z<-optR(A, b, method="cgm")
#' 
#' require(utils)
#' set.seed(129)
#' n <- 7 ; p <- 2
#' X <- matrix(rnorm(n * p), n, p) # no intercept!
#' y <- rnorm(n)
#' Z<-optR(X, y, method="LU")
optR.default<-function(x, y=NULL, method=c("gauss, LU, gaussseidel", "cgm"), iter=500, tol=1e-7, ...){
  
  if(!is.data.frame(x)) x<-data.frame(x)
  if(!is.data.frame(y)) y<-data.frame(y)
  
  # Default fitting
  if(length(method)>1){
    method="LU"
  }
  
  # Assign intial list
  optR<-list()
  
  # Build model
  if(nrow(y)==0){
    # Check for fitting models
    if(method!="LU"){
      warning("b is NULL matrix!!! switching to LU factorization for A decomposition to LU")
    }
    modelf<-as.formula(paste0("~", paste0(colnames(x), collapse="+"), "-1", sep=""))
    optR<-optR(modelf, data=x, method="LU", iter, tol)
  } else
  {
    modelf<-as.formula(paste0(colnames(y), "~", paste0(colnames(x), collapse="+"), "-1", sep=""))
    optR<-optR(modelf, data=cbind.data.frame(x, y), method=method, iter, tol)
  }

  class(optR)<-"optR"
  optR$call<-match.call()
  optR 
}

#' print coefficients for optR class
#' 
#' optR is the default function for optimization
#' @param x : Input of optR class
#' @param ...        : S3 class
#' @export
#' @examples
#' # Solving equation Ax=b
#' A<-matrix(c(6,-4,1, -4,6,-4,1,-4,6), nrow=3,ncol=3, byrow = TRUE)
#' b<-matrix(c(-14,36, 6), nrow=3,ncol=1,byrow=TRUE)
#' Z<-optR(A, b, method="gauss")
#' print(Z)
print.optR<-function(x, ...)
{
  cat("call: \n")
  print(x$call)
  
  # Beta for the coefficients
  if(!is.null(x$beta)){
    cat("\n Coefficients: \n")
    print(x$beta)
  }
}


#' Generate Summary for optR class
#' 
#' summary function generates the summary for the optR class
#' @param object : Input of optR class
#' @param ... : S3 method
#' @export
#' @examples
#' # Solving equation Ax=b
#' A<-matrix(c(6,-4,1, -4,6,-4,1,-4,6), nrow=3,ncol=3, byrow = TRUE)
#' b<-matrix(c(-14,36, 6), nrow=3,ncol=1,byrow=TRUE)
#' Z<-optR(A, b, method="cgm")
#' summary(Z)
summary.optR<-function(object, ...)
{
  # Print the results
  print.optR(object)
  
  # Plot convergence for cgm model
  if(object$method=="cgm"){
    plot(object$conv, xlab="Iterations", ylab="Error")
    lines(object$conv)
    title(main="CGM Convergence Plot")
  }
}