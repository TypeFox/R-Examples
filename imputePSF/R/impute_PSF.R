#' Function to find missing values in time series data and to impute these missing values
#'
#' @param dataIn as inpute time series data with missing values
#' @return returns the time series data with imputed values
#' @export
#' @examples
#' # d <- c(1:5,1:5,1:5,1:5,NA,NA,NA,NA,NA,1:5,1:5,1:5,1:5,NA,NA,NA,NA,NA,1:5,1:5,1:5,1:5)
#' # impute_PSF(d)

#========================================================
# Function "impute_PSF()" starts here-----------------
#========================================================
impute_PSF <- function(dataIn)
{
  a <- missing_patch(dataIn)
  b <- heal_data(dataIn, a)
  return(b)
}

#========================================================
# Function "impute_PSF()" ends here-----------------
#========================================================
