# **************************************************************************************************************
# Functions for permutation
# By default, the number of iterations for bootstrapping, ensembles and permutations is set to 100 
# 
# Functions: 
#      cfPermute:    Main function to call the internal (private) permutation functions and return the results
#      .getPermMatr: Construct the matrix of permuted classes to be used in the permutation iteration
#      .snowRBFperm: Run parallel permutation scripts using the snow/snowfall package in R
# **************************************************************************************************************


cfPermute <- function (inputData, inputClass, bootNum = 100, ensNum = 100, permNum=100, parallel=TRUE, cpus=NULL, type = "SOCK", socketHosts = NULL, progressBar = TRUE, scaling = TRUE) { 
  
  # Get the matrix of permuted classed
  permMatr <- .getPermMatr(inputClass, permNum)
  
  # Start parallelisation of the permutation process
  permObj <- .snowRBFperm(inputData, permMatr, bootNum, ensNum, permNum, parallel, cpus, type, socketHosts, progressBar, scaling)
  
  return (permObj)
}

# Construct the matrix of permuted classes to be used in the permutation iteration
.getPermMatr <- function (classVec, permNum) {
  permMatr = c()
  
  for (permIt in 1:permNum) {
    set.seed(permIt)
    randClass <- as.vector(sample(classVec))
    permMatr  <- rbind(permMatr, randClass)
  }
  
  permMatr <- as.data.frame(permMatr)
  rownames(permMatr) = NULL
  
  return (permMatr)
}


# Run parallel scripts using the snow/snowfall package in R
.snowRBFperm <- function(inputData, permMatr, bootNum, ensNum, permNum, parallel, cpus, type, socketHosts, progressBar, scaling) {
  permList = runTimes = totalTime = avgAcc = c()
  
  # Start the overall timer
  ptm <- proc.time()
  
  tryCatch({
    # Initialisation using given specs from user
    sfInit(parallel=parallel, cpus=cpus, type=type, socketHosts=socketHosts)
    
    # Send the libraries
    sfLibrary("neldermead", character.only=TRUE)
    sfLibrary("e1071",      character.only=TRUE)
    sfLibrary("boot",       character.only=TRUE)
    
    for (k in 1:permNum) { 
      inputClass    <- as.factor(as.matrix(permMatr[k,]))
      
      # Construct the classification ensemble, and count the overall execution time
      execTime      <- system.time(ensRes <- sfLapply(1:ensNum, .boxRadial, inputData, inputClass, bootNum, scaling))
      avgAcc        <- c(avgAcc, round(mean(sapply(ensRes,"[[",1)), digits=2))
      runTimes      <- c(runTimes, execTime[3])
      permList[[k]] <- ensRes
    }
    
    names(runTimes) <- NULL
    totalTime <- proc.time() - ptm
    permList  <- list(avgAcc     = avgAcc,
                      totalTime  = totalTime,
                      execTimes  = round(runTimes, 2),
                      permList   = permList)
    
    class(permList) <- append(class(permList), "cfPermute")
    
    return(permList)
    
    sfStop()
  }, finally=sfStop())
}


.snowRBF <- function(inputData, inputClass, bootNum, ensNum, parallel, cpus, type, socketHosts) {

}
