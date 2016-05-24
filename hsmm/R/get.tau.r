#############################################
# Calculate the length of observations, tau #
#############################################
get.tau <- function(inputData){
  if (length(dim(inputData)) != 2) {
    tau        <- as.integer(length(inputData))
  }
  if (length(dim(inputData)) == 2) {
    tau        <- as.integer(dim(inputData)[2])
  }
  return(tau)
  }
