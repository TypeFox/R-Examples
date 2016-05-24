#' A function to mimic the mode function in Stata.
#' @description This function mimics the functionality of the mode function in 
#' Stata. It does this by calculating the modal category of a vector and 
#' replacing tied categories with a "."to represent a single mode does not exist.
#' @param x a vector, missing values are allowed
#' @param method a character vector of length 1 specifying the way to break ties 
#' in cases where more than one mode exists; either "stata", "sample", or "last". 
#' "stata" provides a "." if more than one mode exists. "sample" randomly samples 
#' from among the tied values for a single mode. "last" takes the final modal 
#' category appearing in the data.
#' @details Specifying method="stata" will result in ties for the mode being 
#' replaced with a "." character. Specifying "sample" will result in the function 
#' randomly sampling among the tied values and picking a single value. Finally, 
#' specifying "last" will result in the function picking the value that appears 
#' last in the original x vector. The default behavior is stata.
#' @return The modal value of a vector if a unique mode exists, else output 
#' determined by method
#' @author Jared E. Knowles
#' @seealso \code{\link{table}} which this function uses
#' @export
#'
#' @examples
#' a <- c(month.name, month.name)
#' statamode(a, method="stata") # returns "." to show no unique mode; useful for ddply
#' statamode(a ,method="sample") # randomly pick one
#' a <- c(LETTERS, "A" , "A")
#' statamode(a)
statamode <- function(x, method = c("last", "stata", "sample")){
  if(is.null(x)){
    return(NA)
  }
  if (missing(method)){
    method <- "stata"
  } else {
    method <- match.arg(method)
  }
  xClass <- class(x)
  x <- as.character(x)
  if (method == 'stata'){
    z <- table(as.vector(x))
    m <- names(z)[z == max(z)]
    if (length(m) == 1){
      if(xClass == "factor"){
        m <- factor(m)
      } else{
        class(m) <- xClass
      }
      return(m)
    }
    return(".")
  }
  else if (method == 'sample'){
    z <- table(as.vector(x))
    m<-names(z)[z == max(z)]
    if (length(m)==1){
      if(xClass == "factor"){
        m <- factor(m)
      } else{
        class(m) <- xClass
      }
      return(m)
    }
    else if (length(m) > 1){
      if(xClass == "factor"){
        m <- factor(m)
      } else{
        class(m) <- xClass
      }
      return(sample(m, 1))
    }
    else if (length(m) < 1){
      if(xClass == "character"){
        return(NA_character_)
      } else if(xClass == "numeric"){
        return(NA_real_)
      } else if(xClass == "integer"){
        return(NA_integer_)
      } else if(xClass == "factor"){
        return(NA_character_)
      }
    }
  }
  else if (method=='last'){
    z <- table(as.vector(x))
    m <- names(z)[z == max(z)]
    if (length(m) == 1){
      if(xClass == "factor"){
        m <- factor(m)
      } else{
        class(m) <- xClass
      }
      return(m)
    }
    else if (length(m) > 1){
      m <- max(x[match(x, m) == max(match(x,m), na.rm=TRUE)], 
               na.rm=TRUE)
      if(xClass == "factor"){
        m <- factor(m)
      } else{
        class(m) <- xClass
      }
      return(tail(m, 1))
    }
    else if (length(m) < 1){
      if(xClass == "character"){
        return(NA_character_)
      } else if(xClass == "numeric"){
        return(NA_real_)
      } else if(xClass == "integer"){
        return(NA_integer_)
      } else if(xClass == "factor"){
        return(NA_character_)
      }
    }
  }
}
