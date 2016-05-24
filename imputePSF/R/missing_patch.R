#' To find out the missing data in given dataset
#'
#' @param dataIn as input time series data with missing values (NAs)
#' @return patch representing table with missing value patches
#' @export

#========================================================
# Function "missing_patch()" starts here-----------------
#========================================================
missing_patch <- function(dataIn)
{
  if(!is.vector(dataIn))
  {
    dataIn <- dataIn[, 1]
  }
  limit1 <- length(dataIn)
  availNA <- which(is.na(dataIn))
  availNA
  m <- diff(availNA)
  m <- append(m,100)
  x <- availNA[1]
  y <- NULL
  for(i in 1:length(availNA))
  {
    if((m[i] != 1) && (!(is.na(m[i] != 1))))
    {
      x <- append(x,availNA[i+1])
      y <- append(y,availNA[i])
    }
  }
  x <- x[!is.na(x)]
  z <- y-x+1
  dataO <- data.frame(x,y,z)
  options( warn = -1 )
  return(dataO)
}

#========================================================
# Function "missing_patch()" ends here-----------------
#========================================================
