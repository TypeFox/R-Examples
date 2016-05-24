################################################################################
##  mvformula: a method to create an object of class mvformula       #
################################################################################

mvformula <- function(x) {
# X: an object that should be coerced to a mvformula
  
	if (is.mvformula(x)) {
		return(x)
    } else if (!inherits(x,"formula"))
    # Turn x into a formula without losing it's environment (eg colnames)
		  x <- formula(x)
    if (is.mvformula(x ))
	return(x )

  resp  <- attr(terms(x),"response")

  if(resp==1){
  # if eval(terms(x)[[2]]) cannot be evaluated, ie if the response is a
  # data.frame this will cause an error
    h <- is.mvabund(eval(terms(x)[[2]]))
  } else
  # if there is no response, the response can't be a mvformula
    h <- FALSE
  if (h)
  {
    class(x) <- c("mvformula", "formula")
    return(x)
  }
  else
  {
  # print warning for objects that might create problems, eg. data.frames
	warning("The response of the formula 'x' is not a mvabund object.
  This can lead to problems in other functions.")  
	   class(x) <- c("mvformula", "formula")
     return(x)
  }
}

################################################################################
# a less stricter version of the function that allows the response in the      #
# formula x to be any kind of matrix                                           #
################################################################################
formulaMva <- function(x){
# X: a formula that is used for mvabund objects
  if ( inherits( try( eval(terms(x)[[2]]) ,silent=TRUE  )[1] , "try-error") |
   inherits( try( model.frame(x) ,silent=TRUE ) [1] , "try-error" ) ) {

    stop("'x' could not be coerced to a 'mvformula'.
		The response and all the independent variables must have the same number of
		cases and the response of the formula must be a matrix or a vector.")

  } else {
      class(x) <- c("mvformula", "formula")
      return(x)
  }
}


################################################################################
# as.mvformula: a function to turn a formula into a formula.mvabund       #
# if the response in x is a data.frame or an unsuitable object the conversion  #
# will fail    !                                                               #
################################################################################
as.mvformula <- function(x)
{  
# X: an object that should be coerced to a mvformula
# note, that the response does not need to be a mvabund object to use this function!
	
  if (is.mvformula(x)) {
		return(x) } else if (!inherits(x,"formula")) {
		object <- as.formula(x) }

  if (is.mvformula(x)) {
		return(x)
  } else if ( inherits( try( eval(terms(x)[[2]]) ,silent=TRUE  )[1] , "try-error")  |
      inherits( try( model.frame(x) ,silent=TRUE )[1], "try-error")  )
      {
        stop("'x' could not be coerced to a 'mvformula'.
		The response and all the independent variables must have the same number of
		cases and the response of the formula must be a matrix or a vector")
      }
      else
      {
          class(x) <- c("mvformula", "formula")
          return(x)
      }
       
}


################################################################################
## a function to see if an object is mvabund	  					 #
################################################################################

is.mvformula <- function(x) {
   inherits(x, "mvformula")
}



# setMethod("show", "mvformula",
# function(object)
# {
# obj<-object@a
# show(obj)
# }
# )    



