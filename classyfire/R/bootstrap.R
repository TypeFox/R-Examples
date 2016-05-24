# ************************************************************************
# Function for bootstrapping
# ************************************************************************

.bootIndx <- function(trainCl, indices) {
  testIndx = trainIndx = c()
  
  # Create the bootTrainClass and bootTestClass based on the passed argument trainClass
  bootTrainCl <- trainCl[indices]
  bootTestCl  <- trainCl[-indices]
  
  # If any of the available classes are not present in bootTrainClass or bootTestClass, 
  # resample by making sure that a sample from each initial class is available in each vector
  
  if (is.element(0, table(bootTrainCl)) || is.element(0, table(bootTestCl))) {   
    classLevels <- levels(trainCl)
    classLength <- nlevels(trainCl)
    
    for (j in 1:classLength) {
      classIndx <- which(trainCl == classLevels[j])
      randIndx  <- sample(classIndx, 2, replace=FALSE)
      
      # Keep at least one sample of each class in bootTrainClass, and similarly in bootTestClass
      testIndx  <- c(testIndx,  randIndx[1])
      trainIndx <- c(trainIndx, randIndx[2])
    }
    
    # Create the bootstrap indices
    indices <- c(trainIndx, sample(setdiff(1:length(trainCl), testIndx), (length(trainCl)-length(trainIndx)), replace=TRUE))
  }
  
  return(indices) 
} 
