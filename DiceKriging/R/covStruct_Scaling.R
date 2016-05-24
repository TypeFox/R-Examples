## --------------
## Scaling class
## --------------

setClass("covScaling",   
         representation(
           d = "integer",           ## (spatial) dimension
           knots = "list",          ## a named list containing the knots
           eta = "list",            ## a list containing the values at knots
           ##knots.n = "integer",   ## number of knots for each dimension
           name = "character",      ## "gauss"
           paramset.n = "integer",  ## number of parameters sets 
           ##   gauss, exp : 1;  powexp : 2
           var.names = "character", ## e.g.  c("Lat", "Long") length d
           ## s.d. of the non-nugget part of error
           sd2 = "numeric",         ## variance (stationarity)
           ## nugget part
           known.covparam = "character",  ## known covariance parameters (except nugget): "All" or "Known"
           nugget.flag = "logical",       ## logical : is there a nugget effect ?
           nugget.estim = "logical",      ## logical : is it estimated (TRUE) or known ?
           nugget = "numeric",            ## nugget (variance)
           ## total number of parameters (except sigma and nugget)
           param.n = "integer"            ## length of knots
         ),
         validity = function(object) {
           
           covset <- c("gauss", "exp", "matern3_2", "matern5_2")
           if (!is.element(object@name, covset)) {
             cat("The list of available covariance functions is:\n", covset, "\n")
             return("invalid character string for 'covtype' argument")
           }
           
           names.knots <- names(object@knots)
           n.knots <- length(object@knots)
           
           if (n.knots>0) {
             if (length(names.knots) == 0) {
               return("the list containing the knots must be named")
             } else if (!all(sapply(object@knots, is.numeric))) {
               return("the knots must be numeric")
             } else if ((!all(names.knots%in%object@var.names)) || (!all(object@var.names%in%names.knots))) {
               return("mismatch between the names of knots and input variables")
             }
             for (i in 1L:n.knots) {
               knots.i <- object@knots[[i]]
               if (any(is.na(knots.i))) {
                 return("missing values not allowed in knots")
               }
               if (is.unsorted(knots.i)) {
                 return("the knots must be sorted")
               }
               if (any(diff(knots.i) <= 0.0)) {
                 return("dupplicated values in knots")
               }
               if (length(knots.i) < 2L) {
                 return("knots must be of length >=2")
               }
             }
             
             n.eta <- length(object@eta)
             if (n.eta>0) {
               if (n.eta != n.knots) {
                 return("the number of values at knots is different from the number of knots")
               } else if (!all(sapply(object@eta, is.numeric))) {
                 return("the values at knots must be numeric")
               } 
               for (i in 1:n.eta) {
                 eta.i <- object@eta[[i]]
                 if (any(eta.i <= 0.0)) {
                   return("the values at knots must be positive")
                 }
               }
             }
           }
           
           if (!identical(object@sd2, numeric(0))) {
             if (object@sd2 < 0) {
               return("The model variance should be non negative")
             }
           }
           
           if (length(object@nugget) > 1) {
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

setMethod("covMatrix", 
          signature = "covScaling", 
          definition = function(object, X, noise.var=NULL) {
            covMatrix(extract.covIso(object), X=scalingFun(X, knots=object@knots, eta=object@eta), noise.var=noise.var)
          }
)

## -----------------------------------------
## Useful METHOD for prediction: covMat1Mat2
## -----------------------------------------

setMethod("covMat1Mat2", 
          signature = "covScaling", 
          definition = function(object, X1, X2, nugget.flag=FALSE) {
            covMat1Mat2(extract.covIso(object), X1=scalingFun(X1, knots=object@knots, eta=object@eta), X2=scalingFun(X2, knots=object@knots, eta=object@eta), nugget.flag=nugget.flag)
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
          signature = "covScaling", 
          definition = function(object){
            param <- unlist(object@eta)
            return(as.numeric(param))
          }
)

setMethod("vect2covparam", 
          signature = "covScaling", 
          definition = function(object, param){
            if (length(param)>0) {
              knots.n <- sapply(object@knots, length)
              ind <- rep(names(knots.n), times=knots.n)
              df <- data.frame(values=param, ind=ind)
              object@eta <- unstack(df)
            }
            return(object)
          }
)

setMethod("coef", 
          signature = signature(object = "covScaling"), 
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
          signature = "covScaling", 
          definition = function(object, X){
            knots.n <- sapply(object@knots, length)
            object.tp <- as(extract.covIso(object), "covTensorProduct")
            bounds <- covParametersBounds(object=object.tp, X=X)
            bounds <- list(lower=rep(1/bounds$upper, times=knots.n), upper=rep(1/bounds$lower, times=knots.n))
            return(bounds)
          }
)

setMethod("paramSample",
          signature = "covScaling", 
          definition = function(object, n, lower, upper, y=NULL, type="all-sd2-nugget"){
            param.n <- object@param.n
            matrixinit <- matrix(runif(n*param.n), nrow = param.n, ncol = n)
            matrixinit <- 1/upper + matrixinit*(1/lower - 1/upper)
            matrixinit <- 1/matrixinit
          }
)



envir.covScaling <- new.env()

setMethod("covMatrixDerivative", 
          signature = "covScaling", 
          definition = function(object, X, C0, k, envir=envir.covScaling) {
            # NOTE : this function MUST be used in a loop over the index k, from 1 to k.max
            if ((k>=1) & (k<=object@param.n)) {
              i <- k
              if (i==1) {
                knots.n <- as.numeric(sapply(object@knots, length))
                k.vec <- rep(1:object@d, times=knots.n)
                l.vec <- sequence(knots.n)
                envir$k.vec <- k.vec
                envir$l.vec <- l.vec
              } else {
                k.vec <- envir$k.vec
                l.vec <- envir$l.vec
              }
              k <- k.vec[i]
              l <- l.vec[i]
              
              if (l==1) {
                object.covTensorProduct <- as(extract.covIso(object), "covTensorProduct")
                fX <- scalingFun(X, knots=object@knots, eta=object@eta)
                Dk <- covMatrixDerivative.dx.covTensorProduct(object=object.covTensorProduct, 
                                                              X=fX, C0=C0, k=k)   
                envir$Dk <- Dk
              } else {
                Dk <- envir$Dk
              }
              df.dkl <- scalingGrad(X=X, knots=object@knots, k)[,l]
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

## ------------
## METHOD show
## ------------

setMethod("show", 
          signature = "covScaling", 
          definition = function(object){
            
            cat("\n")
            cat("Covar. type  :", object@name, ", with scaling \n")
            
            cat("Covar. coeff.")
            if (!identical(object@known.covparam, "All")) cat(", with estimated values for eta")
            cat(":\n")
            
            for (i in 1:object@d) {
              knots.names <- paste("knots", "(", object@var.names[i], ")", sep = "")
              eta.names <- paste("eta", "(", object@var.names[i], ")", sep = "")
              param.names <- c(eta.names, knots.names)
              param.names <- formatC(param.names, width = 12)
              tab <- t(formatC(cbind(object@eta[[i]], object@knots[[i]]), width = 10, digits = 4, format = "f", flag = " "))
              n.i <- length(object@knots[[i]])
              dimnames(tab) <- list(param.names, rep("", n.i))
              print(tab, quote=FALSE)
            }
            
            
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
