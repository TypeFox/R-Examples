
##***********************************************************************
##' S4 class representing a Tensor Sum (TS) covariance kernel.
##'
##' @title S4 class representing a Tensor Sum (TS) covariance kernel
##' @slots d Object of class \code{"integer"} the spatial
##' dimension or number of inputs of the covariance.
##' @docType methods
##' @rdname covTS-methods
##***********************************************************************
setClass("covTS", 		
         representation(
                        d = "integer",                 ## (spatial) dimension
                        inputNames = "character",      ## spatial var names length d
                        kernel = "covMan",             ## "matern5_2add0"
                        kernParNames  = "character",   ## depending on kernel
                        kernParCodes = "integer",      ## 0 fixed, 1 'iso', 2 1par/input 
                        par  = "numeric",              ## params values
                        parN = "integer",              ## number of par
                        parInput = "integer",          ## input number or 0 if 'iso'
                        ## parBlockNames = "character", ## parNames and rep "vars")
                        parLower = "numeric",          ## lower bound on par's
                        parUpper = "numeric",          ## upper bound on par's
                        parBlock = "integer"           ## pointer to param length parN
                        )
)

##*******
# creator
##*******

covTS <- function(inputs = paste("x", 1L:d, sep = ""),
                  d = length(inputs), kernel = "k1matern5_2",
                  dep = NULL, value = NULL, var = 1, ...) {
    
    if (missing(d) & missing(inputs)) stop("at least one of 'd' or 'inputs' must be provided")
    if (length(inputs) != d) stop("'d' must be equal to 'length(inputs)'")
    
    inputNames <- inputs
    
    ## get the kernel as in 'match.fun'. 'kernel' can be a character or
    ## symbol
    if (is.character(kernel) && (length(kernel) == 1L)) {
        if (kernel %in% GPkernNames) kernel <- get(kernel)
        else kernel <- get(kernel, envir = parent.frame(1))
        if (!is(kernel, "covMan")) 
            stop(gettextf("found non-covMan '%s'", kernel), domain = NA)
    } else if (!is(kernel, "covMan")) {
        ## 'kernel' is a symbol???
        if (!is(kernel, "covMan")) {
            if (is.symbol(kernel)) {
                kernel <- eval(substitute(kernel))
            } else {
                stop(gettextf("'%s' is not a covMan object, character or symbol", 
                              deparse(kernel)), domain = NA)
            }
            kernel <- get(as.character(kernel), mode = "any", envir = parent.frame(1))
            if (!is(kernel, "covMan")) 
                stop(gettextf("found non-covMan '%s'", kernel), domain = NA)
        }
    }
    
    if (kernel@d != 1L) {
        stop("'kernel' must specify a one-dimensional kernel",
             " with slot 'd' equal to 1L")
    }
    
    kernName <- kernel@label
    
    ## values defined for the kernel
    DKparNames <- kernel@kernParNames
    
    ## check that "var" is in parNames
    m <- match("var", DKparNames)
    if (is.na(m) || length(m) != 1L) {
        stop("parNames(kernel) must contain one element \"var\"")
    }

    ## remove "var" 
    DKparNames <- DKparNames[-m]
    DKparN <- kernel@parN - 1L
    DKparLower <- kernel@parLower[DKparNames]
    DKparUpper <- kernel@parUpper[DKparNames]
    
    ## caution
    if (DKparN == 1L) names(DKparLower) <-  names(DKparUpper) <- DKparNames
   
    depC <- rep("cst", DKparN)
    names(depC) <- DKparNames
    if (!is.null(dep)) {
        dep <- unlist(dep)
        if ( is.null(names(dep)) || !all(names(dep) %in% DKparNames) ) {
            stop("'dep' must have names corresponding to the kernel parameters")
        }
        depC[names(dep)] <- dep
    } 
    codes <- c("cst" = 1L, "input" = 2L)
    kernParCodes <- rep(1L, DKparN)
    names(kernParCodes) <- DKparNames
    parNames <- character(0)
    parBlock <- integer(0)
    parInput <- integer(0)
    par <- numeric(0)
    parLower <- numeric(0)
    parUpper <- numeric(0)
    i <- 1L
    
    ##======================================================================
    ## prepare information on kernel parameters. This ould be done
    ## without any loop, but would be much less readable
    ## =====================================================================
    for (nm in DKparNames) {
        depi <- codes[depC[nm]]
        if ( nm %in% names(value) ) vali <- value[[nm]]
        else vali <- NA
        if (depi == 1L) {
            if (length(vali) > 1L)
                warning("only the first provided value for ", nm , "will be used.") 
            par <- c(par, vali[1])
            parLower <- c(parLower, DKparLower[nm])
            parUpper <- c(parUpper, DKparUpper[nm])
            parNames <- c(parNames, paste(nm, "_iso_", sep = "."))
            parBlock <- c(parBlock, i)
            parInput <- c(parInput, 0L)
        } else if (depi == 2L) {
            ## if (length(vali) != d)
            ## too much verbosity : removed
            ## warning("Provided vector of values for \"", nm ,
            ##         "\" not of expected length. Recycled.") 
            par <- c(par, rep(vali, length.out = d))
            parLower <- c(parLower, rep(DKparLower[nm], d))
            parUpper <- c(parUpper, rep(DKparUpper[nm], d))
            parNames <- c(parNames, paste(nm, inputNames, sep = "."))
            parBlock <- c(parBlock, rep(i, d))
            parInput <- c(parInput, 1L:d)
        }
        kernParCodes[i] <- depi
        i <- i + 1L
    }
    ## add vars
    if (!missing(var)) {
        var <- rep(var, length.out = d)
    } else var <- rep(1, d)
    ## XXX to be changed?
    wNames <- paste("var", inputNames, sep = ".")
    wLower <- rep(0, d)
    names(wLower) <- wNames
    wUpper <- rep(Inf, d)
    names(wUpper) <- wNames
    par <- c(par, var)
    parLower <- c(parLower, wLower)
    parUpper <- c(parUpper, wUpper)
    parNames <- c(parNames, wNames)
    parBlock <- c(parBlock, rep(DKparN + 1L, d))
    parInput <- c(parInput, 1L:d)
    names(par) <- parNames
    myCov <- new("covTS",
                 d = as.integer(d),
                 inputNames = inputNames,
                 kernel = kernel,
                 kernParNames = DKparNames,
                 kernParCodes = kernParCodes,
                 par = par,
                 parN = length(par),
                 parInput = parInput,
                 parLower = parLower,
                 parUpper = parUpper,
                 parBlock = parBlock)
    myCov
}
##***********************************************************************
## Method to map the parameter vector with the kernel parameters.
##***********************************************************************
setMethod("parMap",
          signature = signature(object = "covTS"),
          definition = function(object, ...){
            
            type <-  c(object@kernParNames, "var")
            kernParCodes <- c(object@kernParCodes, var = 2L)
            nms <- c(object@kernParNames, "var")
            parBlock <- nms[object@parBlock]
            
            valM <- matrix(NA, nrow = object@d, ncol = length(type),
                           dimnames = list(object@inputNames, type))
            d <- object@d
            iPar <- 1L
            for (i in 1L:length(type)) {
              nm <- type[i]
              if (kernParCodes[nm] == 1L) {
                valM[ , nm] <- rep(iPar, d)
                iPar <- iPar + 1L
              } else {
                valM[ , nm] <- iPar - 1L + 1L:d
                iPar <- iPar + d             
              }
            }
            return(valM)
          })
##' npar method for class "TS".
##'
##' npar method for the "covTS" class
##'
##' @param object An object with class "covTS"
##' @return The number of free parmaeters in a covTS covariance.
##' @docType methods
##' @rdname covTS-methods
##'
setMethod("npar",
          signature = signature(object = "covTS"),
          definition = function(object,  ...){
            object@parN
          })


##***********************************************************************
## Utility function to shape some slots containing information about
## parameters. The information can be displayed in list and in matrix
## form as well as in vector form. 
##***********************************************************************
shapeSlot <- function(object, slotName = "par", type = "all", as = "vector"){         
  
    if (identical(type, "all")) {
        if (as == "vector") {
            return(slot(object, slotName))
        }
        ## caution: changing 'type'
        type <-  c(object@kernParNames, "var")
    } else {
        if ( any(duplicated(type)) ||
            !all(type %in% c(object@kernParNames, "var")) ){
            stop("'type' must contain distinct elements, all",
                 "in 'object@kernParNames' or equal to \"var\"")
        }  
    }
    kernParCodes <- c(object@kernParCodes, var = 2L)
    nms <- c(object@kernParNames, "var")
    parBlock <- nms[object@parBlock]
    
    if (as == "vector") {
        ind <- parBlock %in% type
        val <- slot(object, slotName)[ind]
        return(val)
    } else if (as == "list") {
        nms2 <- unlist(lapply(strsplit(names(slot(object, slotName)), "\\."),
                              function(x) paste(x[-1])))
        valL <- list()
        for (i in 1L:length(type)) {
            ind <- parBlock %in% type[i]
            valL[[type[i]]] <- slot(object, slotName)[ind]
            names(valL[[type[i]]]) <- nms2[ind]
        }
        return(valL)
    } else if (as == "matrix") {
        valM <- matrix(NA, nrow = object@d, ncol = length(type),
                       dimnames = list(object@inputNames, type))
        for (i in 1L:length(type)) {
            nm <- type[i]
            ind <- parBlock %in% nm
            val <- slot(object, slotName)[ind]
            if (kernParCodes[nm] == 1L) {
                valM[ , nm] <- rep(val, object@d)
            } else {
                valM[ , nm] <- val
            }
        }
        return(valM)
    }
}

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
          signature = signature(object = "covTS"), 
          definition = function(object){
              cat("Tensor sum covariance kernel\n")
              cat(sprintf("o Dimension 'd' (nb of inputs): %d\n", object@d))
              cat(sprintf("o Kernel (1D): \"%s\" with parameters: %s\n",
                          object@kernel@label,
                          paste(sprintf("\"%s\"", object@kernParNames),
                                collapse = ", ")))
              str <- c("NO (iso)", "YES")
              str <- str[object@kernParCodes]
              str <- paste(object@kernParNames, str, sep = ": ")
              str <- paste(str, collapse = ", ")
              cat(sprintf("o One parameter by input:\n    %s\n", str))
              cat(sprintf("o Number of parameters: %d\n",
                          object@parN))
              
              cat("o Param. values: \n")
              co <- coef(object, type = "all", as = "matrix")
              print(t(co))
          })

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
          signature = signature(object = "covTS"), 
          definition = function(object, type = "all", as = "vector"){         
              shapeSlot(object, slotName = "par", type = type, as = as)
          })

##***********************************************************************
## 'kernelName'
##
## 
##***********************************************************************
setMethod("kernelName",
          signature = signature(object = "covTS"),
          definition = function(object, ...){
              object@kernel@label
          })


##***********************************************************************
## Replacement method
##
## XXX check validity???
##
## NOT WRITTEN YET
##
##***********************************************************************
setMethod("coef<-", 
          signature = signature(object = "covTS", value = "numeric"),
          definition = function(object, type = "all", as = "vector",
              ..., value){
              if (type != "all" || as != "vector") {
                  stop("at the time only implemented values: \"all\" for 'type' ",
                       "and \"vector\" for 'as'")
              }
              if (length(value) != object@parN) {
                  stop(sprintf("'value' must have length %d", object@parN))
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
          signature = signature(object = "covTS"),
          definition = function(object, type = "all", as = "vector"){
              shapeSlot(object, slotName = "parLower", type = type, as = as)            
          })

setMethod("coefLower<-",
          signature = signature(object = "covTS"),
          definition = function(object, type = "all", as = "vector",
              ..., value){
              if (type != "all" || as != "vector") {
                  stop("at the time only implemented values: \"all\" for 'type' ",
                       "and \"vector\" for 'as'")
              }
              if (length(value) != object@parN) {
                  stop(sprintf("'value' must have length %d", npar(object)))
              }
              object@parLower[] <- value
              object   
          })

setMethod("coefUpper", 
          signature = signature(object = "covTS"),
          definition = function(object, type = "all", as = "vector"){
              shapeSlot(object, slotName = "parUpper", type = type, as = as)            
          })

setMethod("coefUpper<-",
          signature = signature(object = "covTS"),
          definition = function(object, type = "all", as = "vector",
              ..., value){
              if (type != "all" || as != "vector") {
                  stop("at the time only implemented values: \"all\" for 'type' ",
                       "and \"vector\" for 'as'")
              }
              if (length(value) != object@parN) {
                  stop(sprintf("'value' must have length %d", npar(object)))
              }
              object@parUpper[] <- value
              object   
          })


## setMethod("compCoefLower", 
##          signature = signature(object = "covTS", X = "matrix"))
## setMethod("compCoefUpper", 
##           signature = signature(object = "covTS", X = "matrix"))

##***********************************************************************
## Methods of the "covMat" familly
##
## Note that the 'covMatrix' method picks the kernel name in the object
## and passes it to the C code after a suitable packing of the
## parameters.
##
## The C_covMatrix uses a 'param' pointer argument and no size 
## indication for it. Make sure that this works correctly even when
## the number of parameters is 0 as in the WhiteNoise case.
##
##***********************************************************************
setMethod("covMat",
          signature = "covTS", 
          definition = function(object, X, Xnew = NULL,
              compGrad = FALSE, checkNames = TRUE,
              index = -1L, ...) {
              
              isXnew <- !is.null(Xnew)
              X <- as.matrix(X)
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
              
              par <- coef(object)
              npar <- length(par)
              parM <- t(parMap(object) - 1L)
              parM[] <- as.integer(parM)
              
              if (compGrad && (index < 1L || index > npar)) {
                  stop("when 'compGrad' is TRUE, 'index' must",
                       " be between 1 and npar(object)")
              }
              compGrad <- as.integer(compGrad)
              index <- as.integer(index) - 1L

              kernFun <- object@kernel@kernel
              
              rho <- new.env()
              environment(kernFun) <- rho

              if (!isXnew) {
                  Cov <- .Call("covMat_covTS", kernFun, t(X), par, parM, compGrad, index, rho,
                               PACKAGE = "kergp")
              } else { 
                  if (compGrad) stop("Gradient computation not implemented when Xnew!=NULL")
                  Cov <- .Call("covMatMat_covTS", kernFun, t(X), t(Xnew), par, parM, 
                               compGrad, index, rho, PACKAGE = "kergp")
              }    
              return(Cov)
          })


##***********************************************************************
## scores method 
## 
## Note that the scores method can only be used when the weight matrix
## is known.
##***********************************************************************

setMethod("scores",
          signature = "covTS", 
          definition = function(object, X, weights, ...) {
              
              X <- as.matrix(X)
              n <- nrow(X)
              d <- ncol(X)
              if (any(is.na(X))) stop("'X' must not contain NA elements") 
              
              par <- coef(object)
              npar <- length(par)
              parM <- t(parMap(object) - 1L)
              parM[] <- as.integer(parM)
              
              kernFun <- object@kernel@kernel
              
              rho <- new.env()
              environment(kernFun) <- rho
              scores <- .Call("scores_covTS", kernFun, t(X), par, parM, weights, rho,
                              PACKAGE = "kergp")
              
          })
