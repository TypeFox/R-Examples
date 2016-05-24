## -----------------------------
## A class for a general kernel
## -----------------------------

setClass("covUser",   
         representation(
           kernel = "function",
           nugget.flag = "logical",    ## logical : is there a nugget effect ?
           nugget = "numeric"
         ), 
         validity = function(object) {
           if (length(object@nugget) > 1L) {
             return("'nugget' must be a single non-negative number. For heteroskedastic noisy observations, use 'noise.var' instead.")
           }           
           if (!identical(object@nugget, numeric(0))) {
             if (object@nugget < 0) {
               return("The nugget effect should be non negative")
             }
           }   
           TRUE
         }
)


## -----------------
## METHOD covMatrix
## -----------------

covMatrixcovUser <- function(object, X, noise.var=NULL) {
  
  d <- ncol(X)
  n <- nrow(X)
  
  C <- matrix(NA, n, n)
  for(i in 1L:n){
    for(j in i:n){
      C[j,i] <- C[i,j] <- object@kernel(X[i,], X[j,])
    }
  }
  
  if (object@nugget.flag) {
    vn <- rep(object@nugget, n)
    C <- C + diag(vn, nrow = n)
  } else if (length(noise.var)>0) {
    vn <- noise.var
    C <- C + diag(noise.var, nrow = n)
  } else {
    vn <- rep(0, n)
  }
  
  return(list(C=C, vn=vn))  
}


setMethod("covMatrix", 
          signature = "covUser", 
          definition = function(object, X, noise.var=NULL) {
            covMatrixcovUser(object=object, X=X, noise.var=noise.var)
          }
)

## -----------------------------------------
## Useful METHOD for prediction: covMat1Mat2
## -----------------------------------------

covMat1Mat2.covUser <- function(object, X1, X2) {
  
  # X1 : matrix n1 x d - containing training points
  # X2 : matrix n2 x d - containing test points
  
  n1 <- nrow(X1)
  n2 <- nrow(X2)
  d <- ncol(X1)
  
  M <- matrix(NA, n1, n2)
  for(i in seq(1,n1)){
    for(j in seq(1,n2)){
      M[i,j] <- object@kernel(X1[i,], X2[j,])
    }
  }
  
  return(M)		
}

setMethod("covMat1Mat2", 
          signature = "covUser", 
          definition = function(object, X1, X2, nugget.flag=FALSE) {
            M <- covMat1Mat2.covUser(object, X1, X2)
            
            if ((!nugget.flag) | (!object@nugget.flag)) {
              return(M)
            } else {
              n1 <- nrow(X1)
              n2 <- nrow(X2)
              out <- .C("C_covMat1Mat2", 
                        as.double(X1), as.integer(n1),
                        as.double(X2), as.integer(n2), 
                        as.integer(ncol(X1)),
                        as.double(0), as.double(object@nugget), "whitenoise",
                        ans = double(n1 * n2), PACKAGE="DiceKriging")
              N <- matrix(out$ans, n1, n2)
              return(M+N) 
            }
          }
)


# ------------------------------
# Useful METHODS for estimation 
# * covparam2vect
# * vect2covparam
# * covParametersBounds
# * covMatrixDerivative
# ------------------------------

setMethod("coef", 
          signature = signature(object = "covUser"), 
          definition = function(object, type = "all", as.list = TRUE){
            val <-
              switch(type,
                     "all" = list(kernel = object@kernel, nugget = object@nugget),
                     "all-nugget"= list(kernel = object@kernel),
                     "kernel" = object@kernel,
                     "nugget" = object@nugget,
                     NULL
              )
            if (as.list) return(val)
            else return(unlist(val, use.names = FALSE))
          }
)


# ------------
# METHOD show
# ------------

setMethod("show", 
          signature = "covUser", 
          definition = function(object){
            cat("\n")
            cat("Covar. type  : user type \n")                    
            print(object@kernel)
            
            if (object@nugget.flag) {
              cat("\nNugget effect :", object@nugget)
              cat("\n\n")
            }
          }
)



