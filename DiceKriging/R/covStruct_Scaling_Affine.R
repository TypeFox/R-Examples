## ---------------------
## Affine scaling class
## ---------------------

setClass("covAffineScaling",   
         representation(
           d = "integer",            ## (spatial) dimension
           knots = "numeric",        ## 2 nodes 
           eta = "matrix",           ## value at nodes, matrix d x 2
           name = "character",       ## "gauss"
           paramset.n = "integer",   ## number of parameters sets 
           ##   gauss, exp : 1;  powexp : 2
           var.names = "character",  ## e.g.  c("Lat", "Long") length d
           ## s.d. of the non-nugget part of error
           sd2 = "numeric",          ## variance (stationarity)
           ## nugget part
           known.covparam = "character",   ## known covariance param (except nugget): "All" or "Known"
           nugget.flag = "logical",        ## logical : is there a nugget effect ?
           nugget.estim = "logical", 	## logical : is it estimated (TRUE) or known ?
           nugget = "numeric",    		## nugget (variance)
           ## total number of parameters (except sigma and nugget)
           param.n = "integer"			## 2*d
         ),
         validity = function(object) {
           
           covset <- c("gauss", "exp", "matern3_2", "matern5_2")
           if (!is.element(object@name, covset)) {
             cat("The list of available covariance functions is:\n", covset, "\n")
             return("invalid character string for 'covtype' argument")
           }
           
           if (!identical(object@sd2, numeric(0))) {
             if (object@sd2 < 0) {
               return("The model variance should be non negative")
             }
           }
           
           if (length(object@nugget)>1) {
             return("Nugget must be a single non-negative number. For heteroskedastic noisy observations, use noise.var instead.")
           }
           
           if (!identical(object@nugget, numeric(0))) {
             if (object@nugget<0) {
               return("The nugget effect should be non negative")
             }
           }
           TRUE
         }
)


## -----------------
## METHOD covMatrix
## -----------------

setMethod("covMatrix", 
          signature = "covAffineScaling", 
          definition = function(object, X, noise.var=NULL) {
            covMatrix(extract.covIso(object), X=affineScalingFun(X, knots=object@knots, eta=object@eta), noise.var=noise.var)
          }
)


## -----------------------------------------
## Useful METHOD for prediction: covMat1Mat2
## -----------------------------------------

setMethod("covMat1Mat2", 
          signature = "covAffineScaling", 
          definition = function(object, X1, X2, nugget.flag=FALSE) {
            covMat1Mat2(extract.covIso(object), X1=affineScalingFun(X1, knots=object@knots, eta=object@eta), X2=affineScalingFun(X2, knots=object@knots, eta=object@eta), nugget.flag=nugget.flag)
          }
)

# ------------------------------
# Useful METHODS for estimation 
# * covparam2vect
# * vect2covparam
# * covParametersBounds
# * covMatrixDerivative
# ------------------------------

setMethod("covparam2vect", 
          signature = "covAffineScaling", 
          definition = function(object){
            param <- matrix(t(object@eta), 2*object@d, 1)
            return(as.numeric(param))
          }
)

setMethod("vect2covparam", 
          signature = "covAffineScaling", 
          definition = function(object, param){
            if (length(param)>0) {
              object@eta <- t(matrix(param, 2, object@d))
            }
            return(object)
          }
)

setMethod("coef", 
          signature = signature(object = "covAffineScaling"), 
          definition = function(object, type = "all", as.list = TRUE){
            val <-
              switch(type,
                     "all" = list(knots = object@knots, eta = object@eta,
                                  sd2 = object@sd2, nugget = object@nugget),
                     "all-nugget"= list(knots = object@knots, eta = object@eta,
                                        sd2 = object@sd2),
                     "all-sd2-nugget" = list(knots = object@knots, eta = object@eta),
                     "knots" = object@knots,
                     "eta" = object@eta,
                     "sd2" = object@sd2,
                     "nugget" = object@nugget,
                     NULL
              )
            if (as.list) return(val)
            else return(unlist(val, use.names = FALSE))
          }
)



setMethod("covParametersBounds", 
          signature = "covAffineScaling", 
          definition = function(object, X){
            object <- as(extract.covIso(object), "covTensorProduct")
            bounds <- covParametersBounds(object=object, X=X)
            bounds <- list(lower=rep(1/bounds$upper, each=2), upper=rep(1/bounds$lower, each=2))
            return(bounds)
          }
)

setMethod("paramSample",
          signature = "covAffineScaling", 
          definition = function(object, n, lower, upper, y=NULL, type="all-sd2-nugget"){
            param.n <- object@param.n
            matrixinit <- matrix(runif(n*param.n), nrow = param.n, ncol = n)
            matrixinit <- 1/upper + matrixinit*(1/lower - 1/upper)
            matrixinit <- 1/matrixinit
          }
)

covMatrixDerivative.dx.covTensorProduct <- function(object, X, C0, k) {
  ## X  : n x d
  n <- nrow(X)
  d <- ncol(X)
  
  ## beware that the index k  starts at 0 in C language 
  
  param <- covparam2vect(object)
  
  out <- .C("C_covMatrixDerivative_dx", 
            as.double(X), as.integer(n), as.integer(d), 
            as.double(param), as.character(object@name),
            as.integer(k), as.double(C0),
            ans = double(n * n),
            PACKAGE="DiceKriging")
  
  return(matrix(out$ans, n, n))
}

envir.covAffineScaling <- new.env()

setMethod("covMatrixDerivative", 
          signature = "covAffineScaling", 
          definition = function(object, X, C0, k, envir=envir.covAffineScaling) {
            ## NOTE : this function MUST be used in a loop over the index k, from 1 to k.max
            if ((k>=1) & (k<=object@param.n)) {
              i <- k
              k <- floor((i-1)/2)+1
              l <- (i-1)%%2+1
              if (l==1) {
                object.covTensorProduct <- as(extract.covIso(object), "covTensorProduct")
                fX <- affineScalingFun(X, knots=object@knots, eta=object@eta)
                Dk <- covMatrixDerivative.dx.covTensorProduct(object=object.covTensorProduct, 
                                                            X=fX, C0=C0, k=k)
                envir$Dk <- Dk
              } else {
                Dk <- envir$Dk
              }
              df.dkl <- affineScalingGrad(X=X, knots=object@knots, k=k)[,l]
              A <- outer(df.dkl, df.dkl, "-")
              dC <- Dk*A
            } else if (k == object@param.n + 1) {
              dC <- C0 / object@sd2
            } else {
              stop("Wrong value of k")
            }
            return(dC)
          }
)


# ------------------------------------
# Useful METHOD for EGO: covVector.dx
# ------------------------------------

# setMethod("covVector.dx", "covAffineScaling",
# function(object, x, X, c) {
#   object.covTensorProduct <- as(extract.covIso(object), "covTensorProduct")
#   fx <- affineScalingFun(matrix(x, nrow=1), knots=object@knots, eta=object@eta)
#   fX <- affineScalingFun(X, knots=object@knots, eta=object@eta)
#   d <- length(x)
#   dkScaling.dx <- rep(Na, d)
#   dk.dx <- covVector.dx.covTensorProduct(object=object.covTensorProduct, x=fx, X=fX, c=c)
#   for (j in 1:d){
#     dkScaling.dx[j] <- dk.dx[j] * affineScalingGrad(X=X, knots=object@knots, k=j)
#   }
#   return(dkScaling.dx)
# }
# )


# ------------
# METHOD show
# ------------

setMethod("show", 
          signature = "covAffineScaling", 
          definition = function(object){
            
            range.names <- paste("eta", "(", object@var.names, ")", sep = "")
            range.names <- formatC(range.names, width = 12)
            tab <- t(formatC(object@eta, width = 10, digits = 4, format = "f", flag = " "))
            if (identical(object@known.covparam, "All")) {
              dimnames(tab) <- list(c("",""), range.names)
            } else {
              dimnames(tab) <- list(c("  Estimate", "  Estimate"), range.names)
            }
            
            cat("\n")
            cat("Covar. type  :", object@name, ", with affine scaling \n")
            
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
