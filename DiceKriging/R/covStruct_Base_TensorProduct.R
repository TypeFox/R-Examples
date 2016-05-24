## --------------------
## Tensor product class
## --------------------

## covTensorProduct : separable (or tensor product) covariances, depending on 1 set of parameters
## Examples - Gaussian, Exponential, Matern with fixed nu=p+1/2, Power-Exponential
##

setClass("covTensorProduct", 		
         representation(
                        d = "integer",            	## (spatial) dimension
                        name = "character",             ## "powexp"
                        paramset.n = "integer",         ## number of parameters sets 
                        ##   gauss, exp : 1;  powexp : 2
                        var.names = "character",  	## e.g.  c("Lat", "Long") length d
			## s.d. of dor the non-nugget part of error
                        sd2 = "numeric",       		## variance (stationarity)
			## nugget part
                        known.covparam = "character",   ## known covariance parameters (except nugget): "All" or "Known"
                        nugget.flag = "logical",  	## logical : is there a nugget effect ?
                        nugget.estim = "logical", 	## logical : is it estimated (TRUE) or known ?
                        nugget = "numeric",    		## nugget (variance)
  			## total number of parameters (except sigma and nugget)
                        param.n = "integer",            ## range.n + shape.n
  			## range part 
                        range.n = "integer",            ## number of distinct range parms
                        range.names = "character",	## their name (usually "theta")
                        range.val = "numeric",          ## their values
  			## shape part, if any 
                        shape.n = "integer",            ## number of distinct shape parms
                        shape.names = "character",      ## their name ("p", "nu", "alpha", etc.)
                        shape.val = "numeric"           ## their values
                        ),
         validity = function(object) {
           
           covset <- c("gauss", "exp", "matern3_2", "matern5_2", "powexp")
           if (!is.element(object@name, covset)) {
             cat("The list of available covariance functions is:\n", covset, "\n")
             return("invalid character string for 'covtype' argument")
           }
           
           if (!identical(object@sd2, numeric(0))) {
             if (object@sd2 < 0) {
               return("The model variance should be non negative")
             }
           }
           
           if (length(object@nugget) > 1L) {
             return("'nugget' must be a single non-negative number. For heteroskedastic noisy observations, use 'noise.var' instead.")
           }
           
           if (!identical(object@nugget, numeric(0))) {
             if (object@nugget < 0) {
               return("The nugget effect should be non negative")
             }
           }
           
           if (!identical(object@range.val, numeric(0))) {
             if (min(object@range.val) < 0) {
               return("The range parameters must have positive values")
             }
             if (length(object@range.val) != object@d) {
               return("Incorrect number of range parameters")
             }
           }
           
           if (!identical(object@shape.val, numeric(0))) {
             if (min(object@shape.val) < 0) {
               return("The shape parameters must have positive values")
             }
             if (length(object@shape.val) != object@d) {
               return("Incorrect number of shape parameters")
             }
             if (identical(object@name, "powexp") && (max(object@shape.val) > 2)) {
               return("The exponents must be <= 2 for a Power-Exponential covariance")
             }
           }
           TRUE
         }
         )


## -----------------
## METHOD covMatrix
## -----------------
covMatrix.covTensorProduct <- function(object, X, noise.var=NULL) {
  
  d <- ncol(X)
  n <- nrow(X)
  
  param <- covparam2vect(object)
  
  out <- .C("C_covMatrix", 
            as.double(X), as.integer(n), as.integer(d), 
            as.double(param), as.double(object@sd2), as.character(object@name), 
            ans = double(n * n),
            PACKAGE="DiceKriging")
  
  C <- matrix(out$ans, n, n)   # covariance matrix when there is no nugget effect
  
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
          signature = "covTensorProduct", 
          definition = function(object, X, noise.var=NULL) {
            covMatrix.covTensorProduct(object=object, X=X, noise.var=noise.var)
          }
)


## -----------------------------------------
## Useful METHOD for prediction: covMat1Mat2
## -----------------------------------------

covMat1Mat2.covTensorProduct <- function(object, X1, X2, nugget.flag=FALSE) {
  
  # X1 : matrix n1 x d - containing training points
  # X2 : matrix n2 x d - containing test points
  
  # X1 <- as.matrix(X1)
  # X2 <- checkNames(X1, X2)
  
  n1 <- nrow(X1)
  n2 <- nrow(X2)
  d <- ncol(X1)
  
  param <- covparam2vect(object)
  
  out <- .C("C_covMat1Mat2", 
            as.double(X1), as.integer(n1),
            as.double(X2), as.integer(n2), 
            as.integer(d),
            as.double(param), as.double(object@sd2), as.character(object@name),
            ans = double(n1 * n2), PACKAGE="DiceKriging")
  
  M <- matrix(out$ans, n1, n2)
  
  if ((!nugget.flag) | (!object@nugget.flag)) {
    return(M)
  } else {
    out <- .C("C_covMat1Mat2", 
              as.double(X1), as.integer(n1),
              as.double(X2), as.integer(n2), 
              as.integer(d),
              as.double(param), as.double(object@nugget), "whitenoise",
              ans = double(n1 * n2), PACKAGE="DiceKriging")
    N <- matrix(out$ans, n1, n2)
    return(M+N)
  }
  
}

setMethod("covMat1Mat2", 
          signature = "covTensorProduct", 
          definition = function(object, X1, X2, nugget.flag=FALSE) {
            covMat1Mat2.covTensorProduct(object=object, X1=X1, X2=X2, nugget.flag=nugget.flag)
          }
)

# ------------------------------
# Useful METHODS for estimation 
# * coef: will replace covparam2vect
# * coef<-: will replace vect2covparam
# * paramSample (optional)
# * covParametersBounds
# * covMatrixDerivative
# ------------------------------

covparam2vect.fun <- function(object){
  if (object@paramset.n==1) {
    param <- object@range.val
  } else param <- c(object@range.val, object@shape.val)
  return(as.numeric(param))
}


setMethod("covparam2vect", 
          signature = "covTensorProduct", 
          definition = covparam2vect.fun
)


vect2covparam.fun <- function(object, param){
  if (length(param)>0) {
    if (object@paramset.n==1) {
      object@range.val <- param
    } else {	
      range.n <- object@range.n
      object@range.val <- param[1:range.n]
      object@shape.val <- param[(range.n+1):length(param)]
    }
  }
  return(object)
}


setMethod("vect2covparam", 
          signature = "covTensorProduct", 
          definition = vect2covparam.fun
)

setMethod("coef", 
          signature = signature(object = "covTensorProduct"), 
          definition = function(object, type = "all", as.list = TRUE){
            val <-
              switch(type,
                     "all" = list(range = object@range.val, shape = object@shape.val,
                                  sd2 = object@sd2, nugget = object@nugget),
                     "all-nugget"= list(range = object@range.val, shape = object@shape.val,
                                        sd2 = object@sd2),
                     "all-sd2-nugget" = list(range = object@range.val, shape = object@shape.val),
                     "range" = object@range.val,
                     "shape" = object@shape.val,
                     "sd2" = object@sd2,
                     "nugget" = object@nugget,
                     NULL
              )
            if (as.list) return(val)
            else return(unlist(val, use.names = FALSE))
          }
)

# POUR EVITER LES CONFLITS AVEC gplab
# setReplaceMethod("coef",
#           signature = signature(object = "covTensorProduct", value = "numeric"),
#           definition = function(object, type = "all", checkValidity = TRUE,
#                                 ..., value) {
#             
#             ## not done during optims
#             if (checkValidity) {
#               expLength <- switch(type, 
#                                   "all" = object@param.n + 1L + as.integer(object@nugget.flag),
#                                   "all-nugget" = object@param.n + 1L,
#                                   "all-sd2-nugget"= object@param.n,
#                                   "range" = object@range.n, 
#                                   "shape" = object@shape.n, 
#                                   "sd2" = 1L,
#                                   "nugget" = 1L)
#               if (length(value) != expLength ) {
#                 stop("'value' must have length ", expLength)
#               }
#             } 
#             if (type == "all") {
#               End <- 0L
#               object@range.val <- value[End + 1L:object@range.n]
#               End <- End + object@range.n
#               if (length(object@shape.n) && object@shape.n > 0L) {
#                 object@range.val <- value[End + 1L:object@shape.n]
#                 End <- End + object@shape.n
#               }
#               object@sd2 <- value[End + 1L]
#               End <- End + 1L
#               if (object@nugget.flag) {
#                 object@nugget.val <- value[End + 1L]
#               }
#             } else if (type == "all-nugget") {
#               End <- 0L
#               object@range.val <- value[End + 1L:object@range.n]
#               End <- End + object@range.n
#               if (length(object@shape.n) && object@shape.n > 0L) {
#                 object@range.val <- value[End + 1L:object@shape.n]
#                 End <- End + object@shape.n
#               }
#               object@sd2 <- value[End + 1L]
#             } else if (type == "all-sd2-nugget") {
#               End <- 0L
#               object@range.val <- value[End + 1L:object@range.n]
#               End <- End + object@range.n
#               if (length(object@shape.n) && object@shape.n > 0L) {
#                 object@range.val <- value[End + 1L:object@shape.n]
#                 End <- End + object@shape.n
#               }
#             } else if (type == "range") {
#               object@range.val <- value
#             } else if (type == "shape") {
#               object@shape.val <- value
#             } else if (type == "sd2") {
#               object@sd2  <- value
#             } else if (type == "nugget") {
#               object@nugget <- value
#             } else stop("bad 'type value")
#             
#             if (checkValidity) validObject(object)
#             return(object) 
#           }
# )

setMethod("covParametersBounds", 
          signature = "covTensorProduct", 
          definition = function(object, X){
            if (object@paramset.n==1) {
              k <- object@range.n
              lower <- rep(1e-10, k)
              upper <- 2 * diff(apply(X, 2, range))
              upper <- as.vector(upper)
            } else if (identical(object@name, "powexp")) {              
              # coef. order : theta_1, ..., theta_d, p_1, ..., p_d
              lower <- rep(1e-10, object@param.n)
              upper <- 2 * diff(apply(X, 2, range))
              k <- object@shape.n
              upper <- as.vector(c(upper, rep(2, k)))
            } else stop("No default values for covariance parameters bounds, the inputs 'lower' and 'upper' are required")
            return(list(lower=lower, upper=upper))
          }
)

setMethod("paramSample",
          signature = "covTensorProduct", 
          definition = function(object, n, lower, upper, y=NULL, type="all-sd2-nugget"){
            param.n <- object@param.n
            matrixinit <- matrix(runif(n * param.n), nrow=param.n, ncol=n)
            matrixinit <- lower + matrixinit*(upper - lower)
          }
)

covMatrixDerivative.covTensorProduct <- function(object, X, C0, k) {
  ## X : n x d
  n <- nrow(X)
  d <- ncol(X)
  
  ##beware that the index k  starts at 0 in C language 
  
  param <- covparam2vect(object)
  
  k <- as.integer(k)
  
  if ((k >=1) & (k <= object@param.n)) {  # derivative with respect to theta_k
    out <- .C("C_covMatrixDerivative", 
              as.double(X), as.integer(n), as.integer(d), 
              as.double(param), as.character(object@name),
              as.integer(k), as.double(C0),
              ans = double(n * n),
              PACKAGE="DiceKriging")
    return(matrix(out$ans, n, n))
  } else if (k==(object@param.n+1)) {     # derivative with respect to sigma^2
    return(C0 / object@sd2)
  } else {
    stop("Wrong value of k")
  }
}

setMethod("covMatrixDerivative",
          signature = "covTensorProduct", 
          definition = function(object, X, C0, k) {
            covMatrixDerivative.covTensorProduct(object=object, X=X, C0=C0, k=k)
          }
)

# ------------------------------------
# Useful METHOD for EGO: covVector.dx
# ------------------------------------

covVector.dx.covTensorProduct <- function(object, x, X, c) {
  
  n <- nrow(X);
  d <- ncol(X);
  
  param <- covparam2vect(object)
  
  out <- .C("C_covVector_dx", 
            as.double(x), as.double(X), 
            as.integer(n), as.integer(d),
            as.double(param), as.character(object@name), 
            as.double(c), 
            ans = double(n * d))
  
  return(matrix(out$ans, n, d))
  
}

setMethod("covVector.dx", 
          signature = "covTensorProduct", 
          definition = function(object, x, X, c) {
            covVector.dx.covTensorProduct(object=object, x=x, X=X, c=c)
          }
)

# ------------
# METHOD show
# ------------

setMethod("show", 
          signature = "covTensorProduct", 
          definition = function(object){
            
            range.names <- paste(object@range.names, "(", object@var.names, ")", sep = "")
            range.names <- formatC(range.names, width = 12)
            val.mat <- matrix(object@range.val, object@d, 1)
            tab <- t(formatC(val.mat, width = 10, digits = 4, format = "f", flag = " "))
            if (identical(object@known.covparam, "All")) {
              dimnames(tab) <- list("", range.names)
            } else {
              dimnames(tab) <- list("  Estimate", range.names)
            }
            
            if (object@paramset.n==2) {
              shape.names <- paste(object@shape.names, "(", object@var.names, ")", sep = "")
              shape.names <- formatC(shape.names, width=12)
              tab <- rbind(tab, shape.names, deparse.level = 0)
              tab <- rbind(tab, formatC(object@shape.val, width = 10, digits = 4, format = "f", flag = " "))
              if (identical(object@known.covparam, "All")) {
                row.names(tab)[3] <- "          "
              } else {
                row.names(tab)[3] <- "  Estimate"
              }
            }
            
            cat("\n")
            cat("Covar. type  :", object@name, "\n")
            
            cat("Covar. coeff.:\n")
            print(t(tab), quote=FALSE)
            cat("\n")
            
            if (identical(object@known.covparam, "All")) {
              cat("Variance:", object@sd2)
            } else {	         
              cat("Variance estimate:", object@sd2)
            }
            cat("\n")
            
            if (object@nugget.flag) {
              if (object@nugget.estim) {
                cat("\nNugget effect estimate:", object@nugget)
              } else cat("\nNugget effect :", object@nugget)
              cat("\n\n")  
            }
          }
)





