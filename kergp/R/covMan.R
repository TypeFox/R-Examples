## see Chambers 'Software for Data Analysis", p 338
setClassUnion("optCharacter", c("character", "NULL"))

setClass("covMan",   	
         representation(
             kernel = "function",
             hasGrad = "logical", 
             acceptMatrix = "logical",      ## TRUE if kernel admit matrices as argu
             label = "character",
             d = "integer",                 ## (spatial) dimension
             inputNames = "optCharacter",   ## spatial var names length d
             parLower = "numeric",          ## lower bound on par's
             parUpper = "numeric",          ## upper bound on par's
             par  = "numeric",              ## params values
             parN = "integer",              ## number of par
             kernParNames  = "character"    ## depending on kernel
         ),
         validity = function(object){
           if (length(object@kernParNames) != object@parN) {
             stop("Incorrect number of parameter names")
           }
         }
)


## -----------------------------
## createur de la classe covMan
## -----------------------------
covMan <- function(kernel,
                   hasGrad = FALSE,
                   acceptMatrix = FALSE,
                   inputs = paste("x", 1:d, sep = ""),
                   d = length(inputs),
                   parNames,
                   par = NULL,
                   label = "covMan", ...) {
    
  if (is.null(par)) par <- as.numeric(rep(NA, length(parNames)))
  if (missing(d) & missing(inputs)) stop("at least one of 'd' or 'inputs' must be provided")
  if (length(inputs) != d) stop("'d' must be equal to 'length(inputs)'")
  
  new("covMan", 
      kernel = as.function(kernel),
      hasGrad = as.logical(hasGrad),
      acceptMatrix = as.logical(acceptMatrix),
      label = as.character(label),
      d = as.integer(d),
      inputNames = as.character(inputs),
      kernParNames = as.character(parNames),
      par = as.numeric(par),
      parN = length(par),    
      ...)
  
}  # TODO : check that the kernel has 2 (eg : brownian) or 3 arguments (parameterised kernel).


setMethod("covMat",
          signature = "covMan", 
          definition = function(object, X, Xnew, compGrad = FALSE, 
                                checkNames = NULL, index = -1L, ...) {

              isXnew <- !is.null(Xnew)
              X <- as.matrix(X)
              
              if (is.null(checkNames)) {
                  checkNames <- TRUE 
                  if (object@d == 1L) {
                      if (ncol(X) == 1L) {
                          checkNames <- FALSE
                      }  
                  }      
              }
              
              if (checkNames) X <- checkX(object, X = X)
              if (any(is.na(X))) stop("'X' must not contain NA elements")

              if (isXnew){
                  Xnew <- as.matrix(Xnew)
                  if (checkNames) Xnew <- checkX(object, X = Xnew)
                  if (ncol(X) != ncol(Xnew)) stop("'X' and 'Xnew' must have the same number of columns")
                  if (any(is.na(Xnew))) stop("'Xnew' must not contain NA elements") 
              } else {
                  Xnew <- X
              }
              
              if (object@acceptMatrix){
                  Cov <- object@kernel(X, Xnew, coef(object))
                  if (!compGrad) attr(Cov, "gradient") <- c()
                  return(Cov)
              }
              
              compGrad <- as.integer(compGrad)
              index <- as.integer(index) - 1L
              par <- coef(object)  
              kernFun <- object@kernel
              rho <- new.env()
              environment(kernFun) <- rho        
              if (!isXnew){
                  Cov <- .Call("covMat_covMan", kernFun, t(X), par, 
                               compGrad, index, rho, PACKAGE = "kergp")
              } else { 
                  if (compGrad) stop("Gradient computation not implemented when Xnew!=NULL")
                  Cov <- .Call("covMatMat_covMan", kernFun, t(X), t(Xnew), par, 
                               compGrad, index, rho, PACKAGE = "kergp")
              }          
              return(Cov) 
              
          })

## setMethod("sd2",
##           signature = signature(object = )

##***********************************************************************
## CAUTION:  when 'type' is a vector and 'as' is "list" or "matrix"
## elements are returned in the order given by 'type'
## which might differ from the standard parameter order.
##
## o 'type' can be "all", or can be a character vector describing a
##          subset of the union U(kernParNaems, "var")
## 
## o 'as'   can be "vector", "list", or "matrix"
##
##***********************************************************************
setMethod("coef", 
          signature = signature(object = "covMan"), 
          definition = function(object){         
            res <- object@par
            names(res) <- object@kernParNames
            res
          })


##***********************************************************************
## Replacement method
##
## XXX check validity???
##
## NOT WRITTEN YET
##
##***********************************************************************

#setReplaceMethod("coef",
#                 signature = signature(object = "covTensorProduct", value = "numeric"),
#                 definition = function(object, type = "all", checkValidity = TRUE,
#                                       ..., value) {
                   

setMethod("coef<-", 
          signature = signature(object = "covMan", value = "numeric"),
          definition = function(object, ..., value){
            if (length(value) != object@parN) {
              stop(sprintf("'value' has length %d but must have length %d",
                           length(value), object@parN))
            }
            object@par[] <- value
            object
          })

##***********************************************************************
## Methods to get/set the parameter bounds?
## One could set bounds by group: range, shape etc.
##
##***********************************************************************
setMethod("coefLower", 
          signature = signature(object = "covMan"),
          definition = function(object){
            object@parLower            
          })

setMethod("coefLower<-",
          signature = signature(object = "covMan"),
          definition = function(object, ..., value){
            if (length(value) != object@parN) {
              stop(sprintf("'value' must have length %d", npar(object)))
            }
            object@parLower[] <- value
            object   
          })

setMethod("coefUpper", 
          signature = signature(object = "covMan"),
          definition = function(object){
            object@parUpper            
          })

setMethod("coefUpper<-",
          signature = signature(object = "covMan"),
          definition = function(object, ..., value){
            if (length(value) != object@parN) {
              stop(sprintf("'value' must have length %d", npar(object)))
            }
            object@parUpper[] <- value
            object   
          })

##***********************************************************************
## coercion method to cleanly extract the kernel slot
##***********************************************************************

setAs("covMan", "function", function(from) from@kernel)


##***********************************************************************
## scores method 
## 
## Note that the scores method can only be used when the weight matrix
## is known.
##***********************************************************************
setMethod("scores",
          signature = "covMan", 
          definition = function(object, X, weights, ...) {

            X <- as.matrix(X)
            n <- nrow(X)
            d <- ncol(X)
            if (any(is.na(X))) stop("'X' must not contain NA elements") 
            
            par <- coef(object)
            npar <- length(par)
            
            ## DO SOME MORE CHECKS HERE

            if (object@acceptMatrix){
                dCov <- attr(object@kernel(X, X, coef(object)), "gradient")
                scores <- rep(NA, npar)          
                for (i in 1:npar){
                    dC <- as.numeric(dCov[[i]][lower.tri(dCov[[i]], diag = TRUE)]) 
                    scores[i] <- sum(weights * dC)
                }
                return(scores)
                ## with array rather than list we would use something like...
                ## 
                ## lt <- lower.tri(matrix(NA, nrow = n , ncol = n), diag = TRUE)
                ## agFun <- function(mat) sum(weights * mat[lt])
                ## scores <- apply(attr(object@kernel(X, X, coef(object)), "gradient"),
                ##                 MARGIN = 1, FUN = agFun)
            }
            
            kernFun <- object@kernel
            rho <- new.env()
            environment(kernFun) <- rho
            scores <- .Call("scores_covMan", kernFun, t(X), par, weights, rho,
                            PACKAGE = "kergp")
            
          })

##***********************************************************************
## The 'show' method must show the kernel name and parameter structure.
## It should also provide information of the parameterisation of the
## structure itself (sharing of the parameters across inputs).
##
##' show method for class "TS"
##' @aliases show,covTS-method
##'
##' @param object XXX
##' @docType methods
##' @rdname covTS-methods
##'
##***********************************************************************
setMethod("show", 
          signature = signature(object = "covMan"), 
          definition = function(object){
            cat("'User' covariance kernel\n")
            argNames <- names(formals(object@kernel))
            cat(paste("o Description:"), object@label, "\n")
            cat(sprintf("o Dimension 'd' (nb of inputs): %d\n", object@d))
#             cat(paste("o Kernel depending on: \"", 
#                       argNames[1], "\", \"", argNames[2], "\"\n", sep=""))
            cat(paste("o Parameters: ",
                      paste(sprintf("\"%s\"", object@kernParNames),
                              collapse = ", "),
                      "\n",sep = ""))
            cat(sprintf("o Number of parameters: %d\n", object@parN))
            if (object@acceptMatrix) {
                cat("o Accept matrix inputs.\n")
            }
            if (object@hasGrad) {
                cat("o Analytical gradient is provided.\n")
            }
            cat("o Param. values: \n")
            co <- as.matrix(coef(object), nrow=1)
            colnames(co) <- "Value"
            print(co)
          })
