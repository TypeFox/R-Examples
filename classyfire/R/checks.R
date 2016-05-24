# ************************************************************************
# Perform initial checks on input arguments
# Throw warning messages in case this step fails
# ************************************************************************


.initCheck <- function (inputData, inputClass, bootNum, ensNum, parallel, scaling) { 
  # Warn and exit if inputData and inputClass are missing or null
  if (missing(inputData) || is.null(inputData)) {
    stop("Argument \"inputData\" is missing (mandatory field).")
  } 
  if (missing(inputClass) || is.null(inputClass)) {
    stop("Argument \"inputClass\" is missing (mandatory field).") 
  } 
  if (('FALSE' %in% sapply(inputData, is.numeric)) == TRUE) {
    stop("Argument \"inputData\" must contain numeric values.")
  }
  
  # Convert the input datasets into the right format:  
  # inputData into a matrix and inputClass into a factor
  inputMatrix <- as.matrix(as.data.frame(inputData))
  inputClass  <- as.factor(as.matrix(inputClass))
    
  classCounts <- table(inputClass) 
  
  if (length(inputClass) != nrow(inputMatrix)) {
    stop("You must provide a class vector equal in length to the number of rows in the data matrix.")
  } else if (any(classCounts < 3)) {
    stop("You need at least three instances of each class in your class vector.")
  } else {
    return (TRUE)
  }
}

.argsCheck <- function (obj, classType) {
  if (missing(obj) || is.null(obj)) {
    stop(paste("Please provide an input argument of type ", classType,".", sep=""))
  }
  if (!(any(class(obj) == classType))) {
    stop(paste("You have to provide an argument of class ", classType,".", sep=""))
  }
}
