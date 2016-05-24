#' Computes the absolute value.
#' 
#' 
#' Helper Function for xpssCompute. R Implementation of the SPSS \code{ABS} Function. 
#'
#' 
#' @usage computeAbs(x) 
#' @param x an atomic numeric or numeric vector or numeric matrix.
#' @return Returns a vector with the absolute values of the input data.
#' @author Bastian Wiessner
#' @seealso \code{\link{abs}}
#' @examples
#' data(fromXPSS)
#' xpssCompute(x = fromXPSS, variables = "V5",fun = "computeAbs")
#' @keywords internal
#' @export
#' 
computeAbs <- function(x){
  #exceptions
  if(!(is.numeric(x))){
    stop("input data has to be numeric")
  }  
  if(ncol(x) == 1){
    out <- base::abs(x)  
  } else{
    stop("computeAbs is limited to one input variable")
  }
  
  return(out)
}

#' Computes the arc-sine.
#'
#'
#' Helper Function for xpssCompute. R Implementation of the SPSS \code{ARSIN} Function.
#' 
#' 
#' @usage computeArsin(x)
#' @param x an atomic numeric or numeric vector or numeric matrix.
#' @return returns a vector with the arc-sine values of the input data.
#' @author Bastian Wiessner
#' @seealso \code{\link{asin}}
#' @examples 
#' xpssCompute(x = fromXPSS, variables = "V5",fun = "computeArsin")
#' @keywords internal
#' @export

computeArsin <- function(x){
  if(!(is.numeric(x))){
    stop("input data has to be numeric")
  } 
  if(ncol(x) == 1){
    out <- base::asin(x)  
  } else{
    stop("computeArsin is limited to one input variable")
  }
  return(out)
}

#' Computes the arc-tan.
#' 
#' 
#' Helper Function for xpssCompute. R Implementation of the SPSS \code{ARTAN} Function.
#'
#'
#' @usage computeArtan(x) 
#' @param x an atomic numeric or numeric vector or numeric matrix.
#' @return returns a numeric, numeric vector or matrix with the arc-tan values of the input data.
#' @author Bastian Wiessner
#' @seealso \code{\link{atan}}
#' @examples 
#' xpssCompute(x = fromXPSS, variables = "V5",fun = "computeArtan")
#' @keywords internal
#' @export

computeArtan <- function(x){
  
  if(!(is.numeric(x))){
    stop("input data has to be numeric")
  } 
  if(ncol(x) == 1){
    out <- base::atan(x)  
  } else{
    stop("computeArtan is limited to one input variable")
  }
  return(out)
}


#' Computes the cosinus.
#' 
#' 
#' Helper Function for xpssCompute. R Implementation of the SPSS \code{COS} Function.
#'
#'
#' @usage computeCos(x) 
#' @param x an atomic numeric or numeric vector or numeric matrix.
#' @return returns a numeric, numeric vector or matrix with the cosinus values of the input data.
#' @author Bastian Wiessner
#' @seealso \code{\link{cos}}
#' @examples 
#' xpssCompute(x = fromXPSS, variables = "V5",fun = "computeCos")
#' @keywords internal
#' @export

computeCos <- function(x){
  
  if(!(is.numeric(x))){
    stop("input data has to be numeric")
  } 
  if(ncol(x) == 1){
    out <- base::cos(x)  
  } else{
    stop("computeCos is limited to one input variable")
  }
  return(out)
}


#' Computes the exponential
#' 
#' 
#' Helper Function for xpssCompute. R Implementation of the SPSS \code{EXP} Function.
#'
#'
#' @usage computeExp(x) 
#' @param x an atomic numeric or numeric vector or numeric matrix.
#' @return returns a numeric, numeric vector or matrix with the expontential values of the input data.
#' @author Bastian Wiessner
#' @seealso \code{\link{exp}}
#' @examples 
#' xpssCompute(x = fromXPSS, variables = "V5",fun = "computeExp")
#' @keywords internal
#' @export

computeExp <- function(x){
  
  if(!(is.numeric(x))){
    stop("input data has to be numeric")
  }
  if(ncol(x) == 1){
    out <- base::exp(x)  
  } else{
    stop("computeExp is limited to one input variable")
  }
  return(out)
}



#' Computes the logarithm base 10
#' 
#' 
#' Helper Function for xpssCompute. R Implementation of the SPSS \code{LG10} Function.
#'
#'
#' @usage computeLg10(x) 
#' @param x an atomic numeric or numeric vector or numeric matrix.
#' @return returns a numeric, numeric vector or matrix with the logarithm base 10 values of the input data.
#' @author Bastian Wiessner
#' @seealso \code{\link{log}}
#' @examples 
#' xpssCompute(x = fromXPSS, variables = "V5",fun = "computeLg10")
#' @keywords internal
#' @export

computeLg10 <- function(x){
  
  if(!(is.numeric(x))){
    stop("input data has to be numeric")
  } 
  if(ncol(x) == 1){
    out <- base::log10(x)  
  } else{
    stop("computeLg10 is limited to one input variable")
  }
  return(out)
}


#' Computes the logarithm naturalis
#' 
#' 
#' Helper Function for xpssCompute. R Implementation of the SPSS \code{LN} Function.
#'
#'
#' @usage computeLn(x) 
#' @param x an atomic numeric or numeric vector or numeric matrix.
#' @return returns a numeric, numeric vector or matrix with the logarithm naturalis values of the input data.
#' @author Bastian Wiessner
#' @seealso \code{\link{log}}
#' @examples 
#' xpssCompute(x = fromXPSS, variables = "V5",fun = "computeLn")
#' @keywords internal
#' @export

computeLn <- function(x){
  
  if(!(is.numeric(x))){
    stop("input data has to be numeric")
  }
  if(ncol(x) == 1){
    out <- base::log(x)  
  } else{
    stop("computeLn is limited to one input variable")
  }
  return(out)
}


#' Computes the logarithm of the gamma function
#' 
#' 
#' Helper Function for xpssCompute. R Implementation of the SPSS \code{LNGAMMA} Function.
#'
#'
#' @usage computeLngamma(x) 
#' @param x an atomic numeric or numeric vector or numeric matrix.
#' @return returns a numeric, numeric vector or matrix with the logarithm of the gamma function values of the input data.
#' @author Bastian Wiessner
#' @seealso \code{\link{lgamma}}
#' @examples 
#' xpssCompute(x = fromXPSS, variables = "V5",fun = "computeLngamma")
#' @keywords internal
#' @export

computeLngamma <- function(x){
  
  if(!(is.numeric(x))){
    stop("input data has to be numeric")
  }
   if(ncol(x) == 1){
    out <- base::lgamma(x)  
  } else{
    stop("computeLngamma is limited to one input variable")
  }
  return(out)
}


#' Returns the remainder of a division.
#' 
#' 
#' Helper Function for xpssCompute. R Implementation of the SPSS \code{MOD} Function.
#'
#'
#' @usage computeMod(x,modulus = NULL)
#' @param x an atomic numeric or numeric vector or numeric matrix.
#' @param modulus atomic numeric \code{x} is divided by.
#' @return returns a numeric vector including the remainder of \code{x} by \code{modulus}.
#' @author Bastian Wiessner
#' @seealso \code{\link{\%\%}}
#' @keywords internal
#' @examples 
#' xpssCompute(x = fromXPSS, variables = "V5",fun = "computeMod", modulus = 2)
#' @export

computeMod <- function(x,modulus = NULL){
  
  if(is.null(modulus)){
    stop("argument modulus is missing")
  } else{
    if(is.numeric(modulus)){
      if(!(is.numeric(x))){
        stop("input data has to be numeric")
      }
      if(ncol(x) == 1){
        out <-  x %% modulus 
      } else{
        stop("computeMod is limited to one input variable")
      }
    }
    else{
      stop("modulus has to be numeric")
    }
  }
  return(out)
}

#' Rounds values.
#' 
#' 
#' Helper Function for xpssCompute. R Implementation of the SPSS \code{RND} Function.
#'
#'
#' @usage computeRnd(x, digits=2)
#' @param x a numeric or complex vector.
#' @param digits atomic numeric argument, indicating the number of decimal places.
#' @return returns the rounded values of \code{x} with the specified number of decimal places.
#' @author Bastian Wiessner
#' @seealso \code{\link{round}}
#' @keywords internal
#' @examples 
#' xpssCompute(x = fromXPSS, variables = "V5",fun = "computeRnd")
#' 
#' xpssCompute(x = fromXPSS, variables = "V5",fun = "computeRnd",digits=1)
#' @export
#' 

computeRnd <- function(x,digits=2){
  if(!(is.numeric(x))){
    stop("input data has to be numeric")
  }
  if(ncol(x) == 1){
    out <- round(x,digits=digits)
  } else{
    stop("computeRnd is limited to one input variable")
  }
  return(out)
}

#' Computes the sine function
#' 
#' 
#' Helper Function for xpssCompute. R Implementation of the SPSS \code{SIN} Function.
#'
#'
#' @usage computeSin(x) 
#' @param x an atomic numeric or numeric vector or numeric matrix.
#' @return returns a numeric, numeric vector or matrix with the sine values of the input data.
#' @author Bastian Wiessner
#' @seealso \code{\link{sin}}
#' @examples 
#' xpssCompute(x = fromXPSS, variables = "V5",fun = "computeSin")
#' @keywords internal
#' @export

computeSin <- function(x){
  
  if(!(is.numeric(x))){
    stop("input data has to be numeric")
  }
  if(ncol(x) == 1){
    out <- base::sin(x)
  } else{
    stop("computeSin is limited to one input variable")
  }
  return(out)
}

#' Computes the square root
#' 
#' 
#' Helper Function for xpssCompute. R Implementation of the SPSS \code{SQRT} Function.
#'
#'
#' @usage computeSqrt(x) 
#' @param x an atomic numeric or numeric vector or numeric matrix.
#' @return returns a numeric, numeric vector or matrix with the square root values of the input data.
#' @author Bastian Wiessner
#' @seealso \code{\link{sqrt}}
#' @examples 
#' xpssCompute(x = fromXPSS, variables = "V5",fun = "computeSqrt")
#' @keywords internal
#' @export

computeSqrt <- function(x){
  
  if(!(is.numeric(x))){
    stop("input data has to be numeric")
  }
  if(ncol(x) == 1){
    out <- base::sqrt(x)
  } else{
    stop("computeSqrt is limited to one input variable")
  }
  return(out)
}
