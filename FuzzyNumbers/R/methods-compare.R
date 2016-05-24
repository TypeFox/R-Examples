#' 
#' @title
#' Possibility of exceedance
#'
#' @description
#' Determines value of possibility of \eqn{e1>=e2}, the result is from
#' range [0,1]. Value 0 indicates no fulfilment of the operator
#' and 1 indicates complete fulfilment.
#' 
#' @template comparison-citation
#' 
#' @param e1 a PiecewiseLinearFuzzyNumber
#' @param e2 a PiecewiseLinearFuzzyNumber
#'
#' @return Returns a value from range [0,1] indicating the possibility of 
#' exceedance of e2 by e1.
#'
#' @exportMethod possibilityExceedance
#' @docType methods
#' @name possibilityExceedance
#' @family comparison-operators
#' @family PiecewiseLinearFuzzyNumber-method
#' @aliases possibilityExceedance,PiecewiseLinearFuzzyNumber,PiecewiseLinearFuzzyNumber-method
#' 
#' @usage
#' \S4method{possibilityExceedance}{PiecewiseLinearFuzzyNumber,PiecewiseLinearFuzzyNumber}(e1, e2)
#' 
#' @templateVar operation possibilityExceedance(a,b)
#' @template comparison-example-exceedance
setGeneric("possibilityExceedance",
           function(e1,e2) standardGeneric("possibilityExceedance"))

setMethod(
  "possibilityExceedance",
  signature(e1 = "PiecewiseLinearFuzzyNumber", e2 = "PiecewiseLinearFuzzyNumber"),
  function (e1, e2)
  {
    fns = fuzzyNumberUnification(e1,e2)
    e1 = fns[[1]]
    e2 = fns[[2]]
    
    if(e1@a4<=e2@a1){
      return(0)
    }
    else if(e1@a2>=e2@a2){
      return(1)
    }
    else{
      
      #alpha cuts from one to zero 
      e1.values = as.numeric(c(e1@a3, e1@knot.right, e1@a4))
      #alpha cuts must be reversed to be from one to zero
      e2.values = as.numeric(rev(c(e2@a1, e2@knot.left, e2@a2)))
      
      above <- (e1.values>=e2.values)
      index = min(which(above == TRUE))
      alphas = rev(c(0,e1@knot.alpha,1))
      
      #substract one from index to find higher alpha cut
      intersection = intersectionY(
                      e1.values[index], alphas[index],
                      e1.values[index-1], alphas[index-1],
                      e2.values[index-1], alphas[index-1],
                      e2.values[index], alphas[index])
          
        return(intersection)
    }
    
  }
)

#' 
#' @title
#' Necessity of exceedance
#'
#' @description
#' Determines value of necessity of \eqn{e1>=e2}, the result is from
#' range [0,1]. Value 0 indicates no fulfilment of the operator
#' and 1 indicates complete fulfilment.
#' 
#' @template comparison-citation 
#' 
#' @param e1 a PiecewiseLinearFuzzyNumber
#' @param e2 a PiecewiseLinearFuzzyNumber
#'
#' @return Returns a value from range [0,1] indicating the necessity of 
#' exceedance of e2 by e1.
#'
#' @exportMethod necessityExceedance
#' @docType methods
#' @name necessityExceedance
#' @family comparison-operators
#' @family PiecewiseLinearFuzzyNumber-method
#' @aliases necessityExceedance,PiecewiseLinearFuzzyNumber,PiecewiseLinearFuzzyNumber-method
#' 
#' @usage
#' \S4method{necessityExceedance}{PiecewiseLinearFuzzyNumber,PiecewiseLinearFuzzyNumber}(e1, e2)
#'
#' @templateVar operation necessityExceedance(a,b) 
#' @template comparison-example-exceedance
setGeneric("necessityExceedance",
           function(e1,e2) standardGeneric("necessityExceedance"))

setMethod(
  "necessityExceedance",
  signature(e1 = "PiecewiseLinearFuzzyNumber", e2 = "PiecewiseLinearFuzzyNumber"),
  function (e1, e2)
  {
    fns = fuzzyNumberUnification(e1,e2)
    e1 = fns[[1]]
    e2 = fns[[2]]

    if(e1@a2<=e2@a1){
      return(0)
    }
    else if(e1@a1>=e2@a2){
      return(1)
    }
    else{
      #here e1 should have inversed membership value for comparison so ve reverse the values
      e1.values = as.numeric(rev(c(e1@a1, e1@knot.left, e1@a2)))
      #alpha cuts from zero to one 
      e2.values = as.numeric(c(e2@a1, e2@knot.left, e2@a2))
      
      above <- (e1.values>=e2.values)
      index = min(which(above == FALSE))
      alphas = c(0,e1@knot.alpha,1)
      
      intersection = intersectionY(
        e1.values[index], alphas[index],
        e1.values[index-1], alphas[index-1],
        e2.values[index], alphas[index],
        e2.values[index-1], alphas[index-1]
      )
      
      return(intersection)      
    }
    
  }
)

#' 
#' @title
#' Possibility of strict exceedance
#'
#' @description
#' Determines value of possibility of \eqn{e1>e2}, the result is from
#' range [0,1]. Value 0 indicates no fulfilment of the operator
#' and 1 indicates complete fulfilment.
#' 
#' @template comparison-citation 
#' 
#' @param e1 a PiecewiseLinearFuzzyNumber
#' @param e2 a PiecewiseLinearFuzzyNumber
#'
#' @return Returns a value from range [0,1] indicating the strict 
#' possibility of exceedance of e2 by e1.
#'
#' @exportMethod possibilityStrictExceedance
#' @docType methods
#' @name possibilityStrictExceedance
#' @family comparison-operators
#' @family PiecewiseLinearFuzzyNumber-method
#' @aliases possibilityStrictExceedance,PiecewiseLinearFuzzyNumber,PiecewiseLinearFuzzyNumber-method
#' 
#' @usage
#' \S4method{possibilityStrictExceedance}{PiecewiseLinearFuzzyNumber,PiecewiseLinearFuzzyNumber}(e1, e2)
#' 
#' @templateVar operation possibilityStrictExceedance(a,b)
#' @template comparison-example-exceedance
setGeneric("possibilityStrictExceedance",
           function(e1,e2) standardGeneric("possibilityStrictExceedance"))

setMethod(
  "possibilityStrictExceedance",
  signature(e1 = "PiecewiseLinearFuzzyNumber", e2 = "PiecewiseLinearFuzzyNumber"),
  function (e1, e2)
  {
    fns = fuzzyNumberUnification(e1,e2)
    e1 = fns[[1]]
    e2 = fns[[2]]
    
    if(e1@a4<e2@a4){
      return(0)
    }
    else if(e1@a2>e2@a4){
      return(1)
    }
    else{
      
      #inverted to fit alpha cut values
      e1.values = as.numeric(rev(c(e1@a3, e1@knot.right, e1@a4)))
      #does not need to be iverted because it should be inversion of alpha value anyway
      e2.values = as.numeric(c(e2@a3, e2@knot.right, e2@a4))
      
      above <- (e1.values>e2.values)
      index = max(which(above == TRUE))
      alphas = c(0,e1@knot.alpha,1)
      
      intersection = intersectionY(
        e1.values[index], alphas[index],
        e1.values[index+1], alphas[index+1],
        e2.values[index], alphas[index],
        e2.values[index+1], alphas[index+1]
      )
      
      return(intersection)
    }
    
  }
)

#' 
#' @title
#' Necessity of strict exceedance
#'
#' @description
#' Determines value of necessity of \eqn{e1>e2}, the result is from
#' range [0,1]. Value 0 indicates no fulfilment of the operator
#' and 1 indicates complete fulfilment.
#' 
#' @template comparison-citation
#' 
#' @param e1 a PiecewiseLinearFuzzyNumber
#' @param e2 a PiecewiseLinearFuzzyNumber
#'
#' @return Returns a value from range [0,1] indicating the strict 
#' necessity of exceedance of e2 by e1.
#'
#' @exportMethod necessityStrictExceedance
#' @docType methods
#' @name necessityStrictExceedance
#' @family comparison-operators
#' @family PiecewiseLinearFuzzyNumber-method
#' @aliases necessityStrictExceedance,PiecewiseLinearFuzzyNumber,PiecewiseLinearFuzzyNumber-method
#' 
#' @usage
#' \S4method{necessityStrictExceedance}{PiecewiseLinearFuzzyNumber,PiecewiseLinearFuzzyNumber}(e1, e2)
#' 
#' @templateVar operation necessityStrictExceedance(a,b) 
#' @template comparison-example-exceedance
setGeneric("necessityStrictExceedance",
           function(e1,e2) standardGeneric("necessityStrictExceedance"))

setMethod(
  "necessityStrictExceedance",
  signature(e1 = "PiecewiseLinearFuzzyNumber", e2 = "PiecewiseLinearFuzzyNumber"),
  function (e1, e2)
  {
    fns = fuzzyNumberUnification(e1,e2)
    e1 = fns[[1]]
    e2 = fns[[2]]
    
    if(e1@a2<e2@a3){
      return(0)
    }
    else if(e1@a1>e2@a4){
      return(1)
    }
    else{
      
      #inverted to fit alpha cut values
      e1.values = as.numeric(rev(c(e1@a1, e1@knot.left, e1@a2)))
      #goes ok
      e2.values = as.numeric(c(e2@a3, e2@knot.right, e2@a4))
     
      above <- (e1.values>e2.values)
      index = max(which(above == TRUE))
      alphas = c(0,e1@knot.alpha,1)
      
      intersection = intersectionY(
        e1.values[index], alphas[index],
        e1.values[index+1], alphas[index+1],
        e2.values[index], alphas[index],
        e2.values[index+1], alphas[index+1]
      )
      
      return(intersection)
      
    }
    
  }
)





#' 
#' @title
#' Possibility of undervaluation
#'
#' @description
#' Determines value of possibility of \eqn{e1<=e2}, the result is from
#' range [0,1]. Value 0 indicates no fulfilment of the operator
#' and 1 indicates complete fulfilment.
#' 
#' @template comparison-citation
#' 
#' @param e1 a PiecewiseLinearFuzzyNumber
#' @param e2 a PiecewiseLinearFuzzyNumber
#'
#' @return Returns a value from range [0,1] indicating the possibility of 
#' exceedance of e2 by e1.
#'
#' @exportMethod possibilityUndervaluation
#' @docType methods
#' @name possibilityUndervaluation
#' @family comparison-operators
#' @family PiecewiseLinearFuzzyNumber-method
#' @aliases possibilityUndervaluation,PiecewiseLinearFuzzyNumber,PiecewiseLinearFuzzyNumber-method
#' 
#' @usage
#' \S4method{possibilityUndervaluation}{PiecewiseLinearFuzzyNumber,PiecewiseLinearFuzzyNumber}(e1, e2)
#' 
#' @templateVar operation possibilityUndervaluation(a,b)
#' @template comparison-example-undervaluation
setGeneric("possibilityUndervaluation",
           function(e1,e2) standardGeneric("possibilityUndervaluation"))

setMethod(
  "possibilityUndervaluation",
  signature(e1 = "PiecewiseLinearFuzzyNumber", e2 = "PiecewiseLinearFuzzyNumber"),
  function (e1, e2)
  {
    fns = fuzzyNumberUnification(e1,e2)
    e1 = fns[[1]]
    e2 = fns[[2]]
    
    if(e1@a1>=e2@a4){
      return(0)
    }
    else if(e1@a2<=e2@a3){
      return(1)
    }
    else{
      
      #alpha cuts from one to zero 
      e1.values = as.numeric(c(e1@a1, e1@knot.left, e1@a2))
      #alpha cuts must be reversed to be from one to zero
      e2.values = as.numeric(rev(c(e2@a3, e2@knot.right, e2@a4)))
      
      above <- (e1.values<=e2.values)
      index = min(which(above == FALSE))
      alphas = c(0,e1@knot.alpha,1)
      
      #substract one from index to find higher alpha cut
      intersection = intersectionY(
        e1.values[index], alphas[index],
        e1.values[index-1], alphas[index-1],
        e2.values[index-1], alphas[index-1],
        e2.values[index], alphas[index])
      
      return(intersection)
    }
    
  }
)

#' 
#' @title
#' Necessity of undervaluation
#'
#' @description
#' Determines value of necessity of \eqn{e1<=e2}, the result is from
#' range [0,1]. Value 0 indicates no fulfilment of the operator
#' and 1 indicates complete fulfilment.
#' 
#' @template comparison-citation 
#' 
#' @param e1 a PiecewiseLinearFuzzyNumber
#' @param e2 a PiecewiseLinearFuzzyNumber
#'
#' @return Returns a value from range [0,1] indicating the necessity of 
#' exceedance of e2 by e1.
#'
#' @exportMethod necessityUndervaluation
#' @docType methods
#' @name necessityUndervaluation
#' @family comparison-operators
#' @family PiecewiseLinearFuzzyNumber-method
#' @aliases necessityUndervaluation,PiecewiseLinearFuzzyNumber,PiecewiseLinearFuzzyNumber-method
#' 
#' @usage
#' \S4method{necessityUndervaluation}{PiecewiseLinearFuzzyNumber,PiecewiseLinearFuzzyNumber}(e1, e2)
#' 
#' @templateVar operation necessityUndervaluation(a,b)
#' @template comparison-example-undervaluation
setGeneric("necessityUndervaluation",
           function(e1,e2) standardGeneric("necessityUndervaluation"))

setMethod(
  "necessityUndervaluation",
  signature(e1 = "PiecewiseLinearFuzzyNumber", e2 = "PiecewiseLinearFuzzyNumber"),
  function (e1, e2)
  {
    fns = fuzzyNumberUnification(e1,e2)
    e1 = fns[[1]]
    e2 = fns[[2]]
    
    if(e1@a2>=e2@a4){
      return(0)
    }
    else if(e1@a4<=e2@a3){
      return(1)
    }
    else{
      
      #alpha cuts from one to zero 
      e1.values = as.numeric(c(e1@a3, e1@knot.right, e1@a4))
      #alpha cuts must be reversed to be from one to zero
      e2.values = as.numeric(rev(c(e2@a3, e2@knot.right, e2@a4)))
      cbind(e1.values,e2.values)
      above <- (e1.values<=e2.values)
      index = min(which(above == FALSE))
      alphas = c(0,e1@knot.alpha,1)
      
      #substract one from index to find higher alpha cut
      intersection = intersectionY(
        e1.values[index], alphas[index],
        e1.values[index-1], alphas[index-1],
        e2.values[index-1], alphas[index-1],
        e2.values[index], alphas[index])
      
      return(intersection)
    }
    
  }
)

#' 
#' @title
#' Possibility of strict undervaluation
#'
#' @description
#' Determines value of possibility of \eqn{e1<e2}, the result is from
#' range [0,1]. Value 0 indicates no fulfilment of the operator
#' and 1 indicates complete fulfilment.
#' 
#' @template comparison-citation 
#' 
#' @param e1 a PiecewiseLinearFuzzyNumber
#' @param e2 a PiecewiseLinearFuzzyNumber
#'
#' @return Returns a value from range [0,1] indicating the necessity of 
#' exceedance of e2 by e1.
#'
#' @exportMethod possibilityStrictUndervaluation
#' @docType methods
#' @name possibilityStrictUndervaluation
#' @family comparison-operators
#' @family PiecewiseLinearFuzzyNumber-method
#' @aliases possibilityStrictUndervaluation,PiecewiseLinearFuzzyNumber,PiecewiseLinearFuzzyNumber-method
#' 
#' @usage
#' \S4method{possibilityStrictUndervaluation}{PiecewiseLinearFuzzyNumber,PiecewiseLinearFuzzyNumber}(e1, e2)
#' 
#' @templateVar operation possibilityStrictUndervaluation(a,b)
#' @template comparison-example-undervaluation
setGeneric("possibilityStrictUndervaluation",
           function(e1,e2) standardGeneric("possibilityStrictUndervaluation"))

setMethod(
  "possibilityStrictUndervaluation",
  signature(e1 = "PiecewiseLinearFuzzyNumber", e2 = "PiecewiseLinearFuzzyNumber"),
  function (e1, e2)
  {
    fns = fuzzyNumberUnification(e1,e2)
    e1 = fns[[1]]
    e2 = fns[[2]]
    
    if(e1@a1>e2@a2){
      return(0)
    }
    else if(e1@a2<e2@a1){
      return(1)
    }
    else{
      
      #alpha cuts from one to zero 
      e1.values = as.numeric(c(e1@a1, e1@knot.left, e1@a2))
      #alpha cuts must be reversed to be from one to zero
      e2.values = as.numeric(rev(c(e2@a1, e2@knot.left, e2@a2)))
      
      above <- (e1.values<=e2.values)
      index = min(which(above == FALSE))
      alphas = c(0,e1@knot.alpha,1)
      
      #substract one from index to find higher alpha cut
      intersection = intersectionY(
        e1.values[index], alphas[index],
        e1.values[index-1], alphas[index-1],
        e2.values[index-1], alphas[index-1],
        e2.values[index], alphas[index])
      
      return(intersection)
    }
    
  }
)

#' 
#' @title
#' Necessity of strict undervaluation
#'
#' @description
#' Determines value of necessity of \eqn{e1<e2}, the result is from
#' range [0,1]. Value 0 indicates no fulfilment of the operator
#' and 1 indicates complete fulfilment.
#' 
#' @template comparison-citation 
#' 
#' @param e1 a PiecewiseLinearFuzzyNumber
#' @param e2 a PiecewiseLinearFuzzyNumber
#'
#' @return Returns a value from range [0,1] indicating the necessity of 
#' exceedance of e2 by e1.
#'
#' @exportMethod necessityStrictUndervaluation
#' @docType methods
#' @name necessityStrictUndervaluation
#' @family comparison-operators
#' @family PiecewiseLinearFuzzyNumber-method
#' @aliases necessityStrictUndervaluation,PiecewiseLinearFuzzyNumber,PiecewiseLinearFuzzyNumber-method
#' 
#' @usage
#' \S4method{necessityStrictUndervaluation}{PiecewiseLinearFuzzyNumber,PiecewiseLinearFuzzyNumber}(e1, e2)
#' 
#' @templateVar operation necessityStrictUndervaluation(a,b)
#' @template comparison-example-undervaluation
setGeneric("necessityStrictUndervaluation",
           function(e1,e2) standardGeneric("necessityStrictUndervaluation"))

setMethod(
  "necessityStrictUndervaluation",
  signature(e1 = "PiecewiseLinearFuzzyNumber", e2 = "PiecewiseLinearFuzzyNumber"),
  function (e1, e2)
  {
    fns = fuzzyNumberUnification(e1,e2)
    e1 = fns[[1]]
    e2 = fns[[2]]
    
    if(e1@a3>e2@a2){
      return(0)
    }
    else if(e1@a4<e2@a1){
      return(1)
    }
    else{
      
      #alpha cuts from one to zero 
      e1.values = as.numeric(c(e1@a3, e1@knot.right, e1@a4))
      #alpha cuts must be reversed to be from one to zero
      e2.values = as.numeric(rev(c(e2@a1, e2@knot.left, e2@a2)))
      cbind(e1.values,e2.values)
      
      above <- (e1.values>e2.values)
      index = max(which(above == FALSE))
      alphas = c(0,e1@knot.alpha,1)
      
      #substract one from index to find higher alpha cut
      intersection = intersectionY(
        e1.values[index], alphas[index],
        e1.values[index+1], alphas[index+1],
        e2.values[index+1], alphas[index+1],
        e2.values[index], alphas[index])
      
      return(intersection)
    }
    
  }
)





#function that unifies number of alpha cuts and their values for any two PiecewiseLinearNumbers
setGeneric("fuzzyNumberUnification",
           function(e1,e2) standardGeneric("fuzzyNumberUnification"))


setMethod(
  "fuzzyNumberUnification",
  signature(e1 = "PiecewiseLinearFuzzyNumber", e2 = "PiecewiseLinearFuzzyNumber"),
  function (e1, e2)
  {
    knot.alpha <- unique(sort(c(e1@knot.alpha, e2@knot.alpha)))
    knot.n <- length(knot.alpha)
    
    if (!isTRUE(all.equal(e1@knot.alpha, knot.alpha)))
      e1 <- as.PiecewiseLinearFuzzyNumber(e1, knot.n=knot.n, knot.alpha=knot.alpha)
    
    if (!isTRUE(all.equal(e2@knot.alpha, knot.alpha)))
      e2 <- as.PiecewiseLinearFuzzyNumber(e2, knot.n=knot.n, knot.alpha=knot.alpha)
          
    return(c(e1,e2))
  }
)

#function for calculating intersection of two lines
setGeneric("intersectionY",
           function(x1,y1,x2,y2,x3,y3,x4,y4) standardGeneric("intersectionY"))

setMethod(
  "intersectionY",
  signature(x1 = "numeric", y1 = "numeric", x2 = "numeric", y2 = "numeric", 
            x3 = "numeric", y3 = "numeric", x4 = "numeric", y4  = "numeric"),
  function (x1,y1,x2,y2,x3,y3,x4,y4)
  {
    x12 = x1 - x2
    x34 = x3 - x4
    y12 = y1 - y2
    y34 = y3 - y4
    
    c = x12 * y34 - y12 * x34
        
    a = x1 * y2 - y1 * x2
    b = x3 * y4 - y3 * x4
    
    #x = (a * x34 - b * x12) / c
    y = (a * y34 - b * y12) / c
        
    return(y)
  }
)
