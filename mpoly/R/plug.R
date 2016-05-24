#' Switch indeterminates in a polynomial
#' 
#' Switch indeterminates in a polynomial
#' 
#' @param p a polynomial
#' @param indeterminate the indeterminate in the polynomial to
#'   switch
#' @param value the value/indeterminate to substitute
#' @return an mpoly object
#' @export
#' @examples
#' 
#' (p <- mp("(x+y)^3"))
#' plug(p, "x", 5)
#' plug(p, "x", "t")
#' plug(p, "x", "y")
#' plug(p, "x", mp("2y"))
#' 
#' plug(p, "x", mp("x + y"))
#' mp("((x+y)+y)^3")
#' 
plug <- function(p, indeterminate, value){
  
  stopifnot(length(value) == 1 || is.mpoly(value))
  stopifnot(length(indeterminate) == 1)
  
  
  ## if plugging in a number
  if(is.numeric(value)){    
    
    pList <- unclass(p)
    pList <- lapply(pList, function(term){
      # term <- pList[[2]]
      indetNdcs <- which(names(term) == indeterminate)
      term[indetNdcs] <- value^unname(term[indetNdcs])
      names(term)[indetNdcs] <- "coef"
      term
    })
    return(mpoly(pList)) # this can be optimized    
    
  }
  
  ## if plugging in a value
  if(is.character(value)) value <- str_trim(value)
  if(is.mpoly(value))     value <- suppressMessages(print.mpoly(value))
  
  if(is.character(value)){
    
    charPoly <- suppressMessages(print.mpoly(p))
    charPoly <- str_replace_all(charPoly, indeterminate, paste0("(",value,")"))
    return(mp(charPoly))# this can be optimized    
    
  }

  
}

