# ************************************************************************
# Function for parallelisation of the classification ensemble 
# ************************************************************************

.snowRBF <- function(inputData, inputClass, bootNum, ensNum, parallel, cpus, type, socketHosts, scaling) {
  tryCatch({
    # Initialisation using given specs from user
    sfInit(parallel=parallel, cpus=cpus, type=type, socketHosts=socketHosts)
    
    # Load the libraries on the workers
    sfLibrary("neldermead", character.only=TRUE)
    sfLibrary("e1071",      character.only=TRUE)
    sfLibrary("boot",       character.only=TRUE)
    sfLibrary("snowfall",   character.only=TRUE)
    
    # Distribute the .boxRadial function to the workers
    # Pass additional arguments: inputData, inputClass, bootNum
    parRBFobj <- sfLapply(1:ensNum, .boxRadial, inputData, inputClass, bootNum, scaling)
    
    # Terminate snowfall
    sfStop()
    
    return(parRBFobj)
  }, finally=sfStop())
}
