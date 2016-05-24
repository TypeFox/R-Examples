## ---------------------------------------------------
## Tensor product class where all parameters are equal
## ---------------------------------------------------


setClass("covIso",   
         representation(
           d = "integer",           ## (spatial) dimension
           name = "character",      ## "gauss"
           paramset.n = "integer",  ## number of parameters sets 
           ##   gauss, exp : 1;  powexp : 2
           var.names = "character", ## e.g.  c("Lat", "Long") length d
           ## s.d. of the non-nugget part of error
           sd2 = "numeric",         ## variance (stationarity)
           ## nugget part
           known.covparam = "character",  ## known covariance parameters (except nugget): "All" or "Known"
           nugget.flag = "logical",    ## logical : is there a nugget effect ?
           nugget.estim = "logical",   ## logical : is it estimated (TRUE) or known ?
           nugget = "numeric",         ## nugget (variance)
           ## total number of parameters (except sigma and nugget)
           param.n = "integer",        ## 1
           ## range part 
           range.names = "character",  ## their name (usually "theta")
           range.val = "numeric"       ## their values
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
           
           if (length(object@nugget) > 1L) {
             return("Nugget must be a single non-negative number. For heteroskedastic noisy observations, use 'noise.var' instead.")
           }
           
           if (!identical(object@nugget, numeric(0))) {
             if (object@nugget < 0) {
               return("The nugget effect should be non negative")
             }
           }
           
           if (!identical(object@range.val, numeric(0))) {
             if (length(object@range.val) > 1L) {
               return("Only one positive value for an isotropic kernel")
             }
             if (min(object@range.val) < 0) {
               return("The range parameter must have positive values")
             }
           }
           TRUE
         }
)


## ----------------
##   FUNCTION - as 
## ----------------

setAs(from = "covIso", to = "covTensorProduct", 
      def = function(from, to) {
        to <- new("covTensorProduct")
        from.names <- slotNames(from)
        for (i in 1:length(from.names)) {
          slot(to, from.names[i]) <- slot(from, from.names[i])
        }
        to@range.n <- from@d
        to@param.n <- to@range.n 
        to@range.val <- rep(from@range.val, from@d)
        return(to)
      }
)


extract.covIso <- function(from) {
  to <- new("covIso")
  from.names <- setdiff(slotNames(from), c("knots", "eta"))
  for (i in 1:length(from.names)) {
    slot(to, from.names[i]) <- slot(from, from.names[i])
  }
  to@range.val <- 1  #rep(1, from@d)
  return(to)
}

## fonction ci-dessous : uniquement pour validation (aucun sens sinon...)
# as.covIso <- function(object){
#   object.iso <- new("covIso")
#   object.names <- setdiff(slotNames(object), c("range.n", "shape.n", "shape.names", "shape.val"))
#   for (i in 1:length(object.names)) {
#     slot(object.iso, object.names[i]) <- slot(object, object.names[i])
#   }
#   object.iso@param.n <- 1
#   object.iso@range.names <- object@range.names[1]
#   object.iso@range.val <- object@range.val[1]
#   return(object.iso)
# }



## -----------------
## METHOD covMatrix
## -----------------

setMethod("covMatrix", 
          signature = "covIso", 
          definition = function(object, X, noise.var=NULL) {
            covMatrix(object=as(object, "covTensorProduct"), X=X, noise.var=noise.var)
          }
)

## -----------------------------------------
## Useful METHOD for prediction: covMat1Mat2
## -----------------------------------------

setMethod("covMat1Mat2", 
          signature = "covIso", 
          definition = function(object, X1, X2, nugget.flag=FALSE) {
            covMat1Mat2(object=as(object, "covTensorProduct"), X1=X1, X2=X2, nugget.flag=nugget.flag)
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
          signature = "covIso", 
          definition = covparam2vect.fun
)

setMethod("vect2covparam", 
          signature = "covIso", 
          definition = vect2covparam.fun
)

setMethod("coef", 
          signature = signature(object = "covIso"), 
          definition = function(object, type = "all", as.list = TRUE){
            val <-
              switch(type,
                     "all" = list(range = object@range.val, sd2 = object@sd2, nugget = object@nugget),
                     "all-nugget"= list(range = object@range.val, sd2 = object@sd2),
                     "all-sd2-nugget" = list(range = object@range.val),
                     "range" = object@range.val,
                     "sd2" = object@sd2,
                     "nugget" = object@nugget,
                     NULL
              )
            if (as.list) return(val)
            else return(unlist(val, use.names = FALSE))
          }
)



setMethod("covParametersBounds",
          signature = "covIso", 
          definition = function(object, X){
            if (object@paramset.n==1) {
              lower <- 1e-10
              upper <- 2 * diff(apply(X, 2, range))
              upper <- max(upper)
            } else stop("No default values for covariance parameters bounds, the inputs 'lower' and 'upper' are required")
            return(list(lower=lower, upper=upper))
          }
)

setMethod("paramSample",
          signature = "covIso", 
          definition = function(object, n, lower, upper, y=NULL, type="all-sd2-nugget"){
            param.n <- object@param.n
            matrixinit <- matrix(runif(n * param.n), nrow=param.n, ncol=n)
            matrixinit <- lower + matrixinit*(upper - lower)
          }
)

setMethod("covMatrixDerivative", 
          signature = "covIso", 
          definition = function(object, X, C0, k) {
            dC <- matrix(0, nrow(C0), ncol(C0))
            object <- as(object, "covTensorProduct")
            if (k == 1) {
              for (j in 1:object@d) {
                dC <- dC + covMatrixDerivative(object=object, X=X, C0=C0, k=j)
              }
            } else if (k==2) {   # derivative with respect to sigma^2
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

setMethod("covVector.dx", 
          signature = "covIso", 
          definition = function(object, x, X, c) {
            covVector.dx.covTensorProduct(object=as(object, "covTensorProduct"), x=x, X=X, c=c)
          }
)

setMethod("show", 
          signature = "covIso", 
          definition = function(object){
            
            range.names <- object@range.names
            range.names <- formatC(range.names, width = 12)
            val.mat <- matrix(object@range.val, 1, 1)
            tab <- t(formatC(val.mat, width = 10, digits = 4, format = "f", flag = " "))
            if (identical(object@known.covparam, "All")) {
              dimnames(tab) <- list("", range.names)
            } else {
              dimnames(tab) <- list("  Estimate", range.names)
            }
            
            cat("\n")
            cat("Covar. type  :", object@name)
            cat(", isotropic")
            cat("\n")
            
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


