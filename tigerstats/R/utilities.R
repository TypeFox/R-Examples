#' Parse formulas
#'
#' utility for extracting portions of formulas.
#'
#' @rdname ParseFormula
#' @keywords internal
#' @usage ParseFormula(formula,...)
#' @param formula, a formula
#' @param \dots additional arguments, should folks decide to add them someday
#' @return an object of class \code{parsedFormula}, used to compute on the language
#' @author Inspired by similar function in package \code{mosaic}.  
#' Included in this package to reduce dependency.
ParseFormula <- function(formula, ...) {
  op <- formula[[1]]
  condition <- NULL
  if (length(formula) == 2) {
    rhs <- formula[[2]]
    lhs <- NULL
  } else if (length(formula) == 3) {
    rhs <- formula[[3]]
    lhs <- formula[[2]]
  } else {
    stop('Invalid formula type.')
  }
  
  if (inherits(rhs, "call") && rhs[[1]] == '|') {
    condition <- rhs[[3]]
    rhs <- rhs[[2]]
  }
  return( structure(list(op=op,lhs=lhs,rhs=rhs,condition=condition), class='parsedFormula') )
}

#' @title Shade Rectangles for Discrete Distributions

#' @description Utility function for pbinomGC ...
#' 
#' @rdname RectShade
#' @keywords internal
#' @usage RectShade(low,high,func,...)
#' @param low lower bound
#' @param high upper bound
#' @param func probability mass function
#' @param \ldots other arguments passed (to modify func)
#' @return graphical side effect only
#' @author Homer White \email{hwhite0@@georgetowncollege.edu}
RectShade <- function(low,high,func,...) { #Utility
  range <- seq(low,high,by=1)
  for (n in range) {
    rect(n-0.5,0,n+0.5,func(n),col="lightblue")
  }
}


#' @title Get a variable from its name

#' @description Primitive utility function, for writing functions that handle formula input.  simpleFind
#' looks first in the environment associated with the data argument.  If nothing is found, it looks
#' in the parent environment, and so on up the chain.  The intent is to allow use of formula constructed
#' from names of variables that may not appear in the data frame of interest, but which are present in the 
#' caller's environment (usually the Global Environment).  Functions that use formulas now are more flexible
#' in an interactive context.
#' 
#' To do:  (1) find a way to make gentler error messages.  (2) as with earlier versions of mosaic, we get
#' the scoping problem when supplied data frame is named 'labels'.
#' 
#' @rdname simpleFind
#' @usage simpleFind(varName,data)
#' @param varName Character string giving the name of the variable to be searched for.
#' @param data Usually a data frame that supplies the some or all of the variables for a formula 
#' that is has been passed to the calling function.
#' @keywords internal
#' @author Homer White \email{hwhite0@@georgetowncollege.edu}
simpleFind <- function(varName,data) {
  
  if (is.null(data)) {
    return(get(varName,inherits=T))
  }
  
  tryCatch({get(varName,envir=as.environment(data))},
           error=function(e) {
             # is data name on the search path?
             dataName <- deparse(substitute(data))
             
             # will throw the error if data is not on search path 
             get(dataName,inherits=T)  
             
             # otherwise, user probably intends that this particular variable
             # is outside the stated data, so look for it:
             possibleVar <- get(varName,inherits=T)
             
             # the following is a bit of a hack, but here goes:
             # sometimes the name coincides with one of R's base functions.
             # Perhaps then the user has the variable in the Globabl Environment.
             # Look for it there:
             if (is.function(possibleVar)) {
               possibleVar <- get(varName,inherits=T,envir = globalenv())
             }
             possibleVar
           }
  )
  
}


#' @title Shade Under Density Curves

#' @description Utility function for ptGC, pnormGC, pchisqGC, possibly others
#' @keywords internal
#' @rdname UnderShade
#' @usage UnderShade(low,high,func,...)
#' @param low lower bound
#' @param high upper bound
#' @param func density function
#' @param \ldots other arguments passed (to modify func)
#' @return graphical side effect only
#' @author Homer White \email{hwhite0@@georgetowncollege.edu}
UnderShade <- function(low,high,func,...) { #Utility
  x.coords <- c(low,seq(low,high,length.out=301),high)
  y.coords <- c(0,func(seq(low,high,length.out=301),...),0)
  polygon(x.coords,y.coords,col="lightblue",cex=2)
}




