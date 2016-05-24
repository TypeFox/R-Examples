
# this function is internal, not meant for the user

# R code to process user input, then call
# the corresponding function of the C++ object


.sl <- function( z, theta0, alpha0 =NULL, method ="clr",
              tol =0.001, output ="T" )  {
# overload by if-else statements
  if( !is.null(alpha0) && !is.numeric(alpha0) )  {
    if( is.numeric(method) ) {
      if( !missing(tol) )  {
        if( is.character(tol) )  output <- tol
          else  stop( "'output' must be \"T\", \"V\" or \"B\"" )
      }
      tol <- method
    }
    method <- alpha0
    alpha0 <- NULL
  }
  method <- toupper(method)
  output <- toupper(output)
  met <- integer(1)
  if( method=="CLR" )  met <- 1  else  {
    if( method=="AF" )  met <- 2  else  {
      if( method=="MC" )  met <- 3  else
        stop( "'method' must be \"CLR\", \"AF\" or \"MC\"" )
    }
  }
  value <- verbose <- logical(1)
  if( output=="T" )  {
    value <- FALSE
    verbose <- TRUE  
  } else {
    if( output=="V" )  {
      value <- TRUE
      verbose <- FALSE
    } else {
      if( output=="B" )
        value <- verbose <- TRUE
      else  stop( "'output' must be \"T\", \"V\" or \"B\"" )
    }
  }
  if(value) {
    result <- double(1)
    result <- if( is.null(alpha0) )
        (z$CppObj)$sl5( met, as.integer(verbose),
          as.integer(value), tol, theta0 )
      else
        (z$CppObj)$sl6( met, as.integer(verbose),
          as.integer(value), tol, theta0, alpha0 )
    return( result )
  }  else  {
    if( is.null(alpha0) | !is.numeric(alpha0) )
      (z$CppObj)$sl3( met, tol, theta0 )
    else
      (z$CppObj)$sl4( met, tol, theta0, alpha0 )
  }
}

