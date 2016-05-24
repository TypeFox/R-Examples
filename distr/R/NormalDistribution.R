
################################
##
## Class: NormParameter
##                             
################################



## Access methods
setMethod("mean", "NormParameter", function(x, ...) x@mean)
setMethod("sd", signature(x = "NormParameter"), function(x, ...) x@sd)
## Replace Methoden
setReplaceMethod("mean", "NormParameter", 
                  function(object, value){ object@mean <- value; object})
setReplaceMethod("sd", "NormParameter", 
                  function(object, value)
                      { object@sd <- as.matrix(value); object})


validNormParameter <- function(object){
    if(!is.matrix(sd(object)) && is.numeric(mean(object))) return(TRUE)
    if(nrow(sd(object)) != ncol((sd(object))))
        stop("Covariance matrix not sqared")
    if(nrow(sd(object)) != length(mean(object)))
        stop("Covariance matrix and mean vector do not have the same dimension")
    return(TRUE)
}

setValidity("NormParameter", validNormParameter)


################################
##
## Class: UniNormParameter
##
################################

setClass("UniNormParameter", contains = "NormParameter")


setValidity("UniNormParameter", function(object){
  if(length(mean(object)) != 1)
    stop("mean has to be a numeric of length 1")    
  if(length(sd(object)) != 1)
    stop("sd has to be a numeric of length 1")    
  sd <- as.numeric(sd(object))
  if(sd <= 0)
    stop("sd has to be positive")
  else return(TRUE)
})


################################
##
## Class: normal distribution
##
################################

Norm <- function(mean = 0, sd = 1) {
   N <- new("Norm", mean = mean, sd = sd)
   N@Symmetry <- SphericalSymmetry(mean)
   N}

## wrapped access methods
setMethod("mean", "Norm", function(x, ...) mean(param(x)))
setMethod("sd", signature(x = "Norm"), function(x) sd(param(x)))
## wrapped replace methods 
setMethod("mean<-", "Norm", 
           function(object, value) new("Norm", mean = value, sd = sd(object)))
setMethod("sd<-", "Norm", 
           function(object, value) new("Norm", mean = mean(object), sd = value))

## clipped moments for normal distribution: found in distrEx...

###setMethod("m1df", "Norm", 
###          function(object){
###            function(t) -d(object)(t) * sd(param(object))^2 
###                        + mean(param(object)) * p(object)(t)
###          })
###
###setMethod("m2df", "Norm", 
###          function(object){
###            mean <- mean(param(object))
###            sd <- sd(param(object))
###            d <- d(object)
###            p <- p(object)
###            function(t) -(t-mean) * d(t) * sd^2 + p(t) * sd^2 - 
###                        2 * mean * d(t) * sd^2 + mean^2 * p(t) 
###          })
###
## convolution operator for normal distributions

setMethod("+", c("Norm","Norm"),
          function(e1,e2){
              N<- new("Norm", sd = sqrt(sd(e1)^2 + sd(e2)^2), 
                   mean = mean(e1) + mean(e2),  .withArith = TRUE)
              N@Symmetry <- SphericalSymmetry(mean(e1)+mean(e2))
              N 
          })

## extra methods for normal distribution
setMethod("+", c("Norm","numeric"),
          function(e1, e2){
            if (length(e2)>1) stop("length of operator must be 1")
            N <- new("Norm", mean = mean(e1) + e2, sd = sd(e1), .withArith = TRUE) 
            N@Symmetry <- SphericalSymmetry(mean(e1)+ e2)
            N           
          })

setMethod("*", c("Norm","numeric"),
          function(e1, e2){
            if (length(e2)>1) stop("length of operator must be 1")
            if (isTRUE(all.equal(e2,0))) 
                return(new("Dirac", location = 0, .withArith = TRUE))
            N<- new("Norm", mean = mean(e1) * e2, 
                 sd = sd(e1) * abs(e2), .withArith = TRUE)
            N@Symmetry <- SphericalSymmetry(mean(e1)*e2)
            N           
          })
