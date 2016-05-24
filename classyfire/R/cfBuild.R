# ************************************************************************
# Main function for the construction of the classification ensemble
# ************************************************************************

cfBuild <- function (inputData, inputClass, ...) UseMethod("cfBuild") 

cfBuild.default <- function(inputData, inputClass, bootNum = 100, ensNum = 100, parallel = TRUE, cpus = NULL, type = "SOCK", socketHosts = NULL, scaling = TRUE, ...) {
  
  if(.initCheck(inputData, inputClass, bootNum, ensNum, parallel, scaling)) {
    
    # Convert the input arguments into the right format
    inputData   <- as.matrix(as.data.frame(inputData))
    inputClass  <- as.factor(as.matrix(inputClass))
    
    # Construct the classification ensemble, and count the overall execution time
    execTime    <- system.time(svmObj <- .snowRBF(inputData, inputClass, bootNum, ensNum, parallel, cpus, type, socketHosts, scaling))
    ensList     <- list(testAcc     = round(sapply(svmObj,"[[", 1), 2),
                        trainAcc    = round(sapply(svmObj,"[[", 2), 2),
                        optGamma    = sapply(svmObj,"[[", 3), 
                        optCost     = sapply(svmObj,"[[", 4),
                        totalTime   = execTime,
                        runTime     = t(sapply(svmObj,"[[", 5))[,3],
                        confMatr    = lapply(svmObj,"[[", 6),
                        predClasses = lapply(svmObj,"[[", 7),
                        testClasses = lapply(svmObj,"[[", 8),
                        missNames   = lapply(svmObj,"[[", 9),
                        accNames    = lapply(svmObj,"[[", 10), 
                        testIndx    = lapply(svmObj,"[[", 11), 
                        svmModel    = lapply(svmObj,"[[", 12))
    
    class(ensList) <- append(class(ensList), "cfBuild")
    
    return(ensList)
  }
}
