#' Transform values of a vector
#' 
#' @author Bruno Vilela
#' 
#' @description Transform each element of a vector.
#' 
#' @param x A vector to be transformed.
#' @param y levels to be transformed. 
#' @param z The value to be atributed to each level (same order as y).
#' @param NUMERIC logical, if \code{TRUE} z will be considered numbers.
#' 
#' @return Return a vector with changed values.
#' 
#' @examples \dontrun{
#' status <- sample(c("EN","VU", "NT", "CR", "DD", "LC"), 30, replace=TRUE) 
#' TE <- "Threatened"
#' NT <- "Non-Threatened"
#' new <- c(TE, TE, NT, TE, "Data Deficient", NT)
#' old <- c("EN","VU", "NT", "CR", "DD", "LC")
#' statustrans <- lets.transf(status, old, new, NUMERIC=FALSE)
#' 
#' }
#' 
#' @export



lets.transf <- function (x, y, z, NUMERIC = TRUE) 
{
  if(is.factor(x)){
    x <- as.numeric(levels(x))[x]
  }
  
  if(is.factor(y)){
    y <- as.numeric(levels(y))[y]
  }
  
  for(i in 1:length(y)){  
    x[x == y[i]] <- z[i]
  }
  if(NUMERIC){
    x <- as.numeric(x)
  }
  return(x)
  
}
