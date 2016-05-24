#******************************************************************************#
# Verify and pre-process inputs                                                #
#******************************************************************************#
#                                                                              #
# Inputs                                                                       #
#                                                                              #
#  data.x         an object of class data.frame.                               #
#                 The structure of the data.frame must be                      #
#                 \{patient ID, event time, event indicator\}.                 #
#                 Patient IDs must be of class integer or be able to be        #
#                 coerced to class integer without loss of information.        #
#                 Missing values must be indicated as NA.                      #
#                                                                              #
#  data.z         an object of class data.frame.                               #
#                 The structure of the data.frame must be                      #
#                 \{patient ID, time of measurement, measurement(s)\}.         #
#                 Patient IDs must be of class integer or be able to be        #
#                 coerced to class integer without loss of information.        #
#                 Missing values must be indicated as NA.                      #
#                                                                              #
# Outputs                                                                      #
#                                                                              #
#  Return a list                                                               #
#                                                                              #
#  data.x        Same as input with: ids coerced to integer; NAs removed;      #
#                                                                              #
#  data.z        Same as input with: ids coerced to integer; NAs removed;      #
#                missing data cases set to 0.                                  #
#                                                                              #
#******************************************************************************#
preprocessInputs <- function(data.x, data.z) {

  #--------------------------------------------------------------------------#
  # Verify sufficient number of columns in datasets                          #
  #--------------------------------------------------------------------------#
  nc <- ncol(data.x)
  if( nc != 3L ) stop("data.x must include {ID, time, delta}.")

  ncz <- ncol(data.z)
  if( ncz < 3L ) stop("data.z must include {ID, time, measurement}.")

  #--------------------------------------------------------------------------#
  # ensure that patient ids are integers                                     #
  #--------------------------------------------------------------------------#
  if( !is.integer(data.z[,1L]) ) {
    data.z[,1L] <- as.integer(round(data.z[,1L],0))
    cat("Patient IDs in data.z were coerced to integer.\n")
  }

  if( !is.integer(data.x[,1L]) ) {
    data.x[,1L] <- as.integer(round(data.x[,1L],0))
    cat("Patient IDs in data.x were coerced to integer.\n")
  }

  #--------------------------------------------------------------------------#
  # Remove any cases for which all covariates are NA                         #
  #--------------------------------------------------------------------------#
  rmRow <- apply(data.z, 1, function(x){all(is.na(x))})
  data.z <- data.z[!rmRow,]

  #--------------------------------------------------------------------------#
  # Set missing cases to 0.0                                                 #
  #--------------------------------------------------------------------------#
  tst <- is.na(data.z)
  data.z[tst] <- 0.0

  #--------------------------------------------------------------------------#
  # Remove any cases for which response is NA                                #
  #--------------------------------------------------------------------------#
  tst <- is.na(data.x[,2L])
  data.x <- data.x[!tst,]

  #--------------------------------------------------------------------------#
  # Determine if range of data.z is (0,1)                                    #
  #--------------------------------------------------------------------------#
  if( any(data.z[,2L] < {-1.5e-8}) ) {
    stop("Time is negative in data.z.")
  }

  if( any(data.x[,2L] < {-1.5e-8}) ) {
    stop("Time is negative in data.x.")
  }

  return(list(data.x = data.x,
              data.z = data.z))

}
