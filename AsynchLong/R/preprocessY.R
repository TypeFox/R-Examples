preprocessY <- function(data.y) {

  #------------------------------------------------------------------#
  # Verify sufficient number of columns in datasets                  #
  #------------------------------------------------------------------#
  ncy <- ncol(data.y)
  if( ncy < 3L ) {
    stop("data.y must include {ID, time, measurement}.",
         call. = FALSE)
  }

  #------------------------------------------------------------------#
  # ensure that patient ids are integers                             #
  #------------------------------------------------------------------#
  if( !is.integer(data.y[,1L]) ) {
    data.y[,1L] <- as.integer(round(data.y[,1L],0))
    cat("Patient IDs in data.y were coerced to integer.\n")
  }

  #------------------------------------------------------------------#
  # Remove any cases for which response is NA                        #
  #------------------------------------------------------------------#
  tst <- is.na(data.y[,3L])
  data.y <- data.y[!tst,]

  #------------------------------------------------------------------#
  # Determine if any times in data.y are < 0                         #
  #------------------------------------------------------------------#
  if( any(data.y[,2L] < {-1.5e-8}) ) {
    stop("Time is negative in data.y.", call. = FALSE)
  }

  return(data.y)

}
