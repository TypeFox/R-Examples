#' Load the Input-Output and Final demand tables
#' 
#' This function loads the demand tables
#' and defines all variables for the decomposition
#' 
#' @param x the intermediate demand table, it has dimensions GN x GN (G = no. of country, N = no. of industries),
#'  excluding the first row and the first column which contains the country names,
#'  and the second row and second column which contain the industry names for each country.
#'  In addition, an extra row at the end should contain final demand.
#' @param y the final demand table it has dimensions GN x MN,
#'  excluding the first row and the first column which contains the country names,
#'  the second column which contains the industry names for each country,
#'  and second row which contains the five decomposed final demands (M).
#' @return a decompr class object
#' @author Bastiaan Quast
#' @details Adapted from code by Fei Wang.
#' @export


load_tables <- function(x, y) {
  
  # Part 1: getting the rownames etc.
  secreg <- as.character(x[2,-c(1,2)])
  GN     <- length(secreg)
  secnam <- unique(secreg)
  N      <- length(secnam)
  G      <- as.integer(GN / N)
  regnam <- unique( as.character(x[1,-c(1,2)] ) )
  
  x <- x[ -c(1,2),-c(1,2) ]
  x <- apply( x, 2, as.numeric )
  output <- x[ dim(x)[1], ]
  
  x <- x[ 1:GN, ]
  
  y <- data.matrix( y[3:(GN+2),3:((5*G)+2)] )
  
  warning("The API for the decomp function has changed,
  it now uses load_tables_vectors instead of load_tables,
  for more info see http://qua.st/decompr/decompr-v2/.")
  
  load_tables_vectors(x,
                      y,
                      regnam,
                      secnam,
                      output)
  
}