#' Impute missing values in a variable
#' 
#'  Takes in a vector and a value, and returns the vector with missing values imputed 
#'  with that value
#'  @param vector a vector with missing values
#'  @param value the value to be used for imputation
#'  @details
#'  \code{imputemiss} imputes the missing (NA) values in the vector with a specified value.
#'  The function simplifies the code for imputation.
#'  @return \code{vector} of the same class as input vector with imputed missing values
#'  @author Akash Jain
#'  @seealso \code{\link{decile}}, \code{\link{pentile}}, \code{\link{outliers}}
#'  @examples
#'  # Scores vector
#' scores <- c(1, 2, 3, NA, 4, NA)
#'
#' # Imputd scores vector
#' scoresImp <- imputemiss(vector = scores, value = 5)
#'  @export
imputemiss <- function(vector, value) {
  numNA <- sum(is.na(vector))
  if(numNA == 0) {
    stop('There are no missing values in the vector')
  } else{
    vector[is.na(vector)] <- value
  }
  return(vector)
}