#' Pretty printing of a list of multivariate polynomials.
#' 
#' This function iterates print.mpoly on an object of class
#' mpolyList.
#' 
#' @param x an object of class mpolyList
#' @param varorder the order of the variables
#' @param order a total order used to order the monomials in the
#'   printing
#' @param ... arguments to pass to print.mpoly
#' @usage \method{print}{mpolyList}(x, varorder = vars(x), order,
#'   ...)
#' @return Invisible character vector of the printed objects.
#' @export
#' @examples
#' 
#' mL <- mp(c('x + 1', 'y - 1', 'x y^2 z  +  x^2 z^2  +  z^2  +  x^3'))
#' mL
#' print(mL, order = 'lex')
#' print(mL, order = 'glex')
#' print(mL, order = 'grlex')
#' print(mL, order = 'glex', varorder = c('z','y','x'))
#' print(mL, order = 'grlex', varorder = c('z','y','x'))
#' print(mL, varorder = c('z','y','x'))
#' s <- print(mL, varorder = c('z','y','x'))
#' str(s)
#' 
print.mpolyList <- function(x, varorder = vars(x), order, ...){
  
  stopifnot(is.mpolyList(x))
  n <- length(x)
  vars <- vars(x)
  
  if(!missing(varorder) && !all(vars %in% varorder)){
    error <- paste("if specified, varorder must contain all computed vars - ", 
        paste(vars, collapse = ", "), sep = "")
    stop(error, call. = FALSE)
  }
  
  if(missing(varorder) && !missing(order)){
    message <- paste(
      'using variable ordering - ',
      paste(vars, collapse = ', '),
      sep = ''
    )
    message(message)
  }
  
  if(missing(order)){
    polys <- sapply(x, print, varorder = varorder, ...)  
  } else {
    polys <- sapply(x, print, varorder = varorder, order = order, ...)    
  }
 
  invisible(polys)
}








