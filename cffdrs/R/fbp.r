fbp  <- function(input = NULL, output = "Primary", m = NULL, cores = 1){  
  #############################################################################
  # Description:
  #   An internal function used to setup the calculation of the Fire Behavior 
  #   Prediction (FBP) system.This function moves the logic of the FBP system
  #   equations into FBPcalc.R and sets up the use of that function here.
  #
  #
  # Args:
  #   input:  Data frame of required and optional information needed to 
  #           calculate FBP function. View the arguments section of the fbp 
  #           manual (fbp.Rd) under "input" for the full listing of the 
  #           required and optional inputs.
  #   output: What fbp outputs to return to the user. Options are "Primary", 
  #           "Secondary" and "All".
  #   m:      Optimal number of pixels at each iteration of computation.
  #   cores:  Number of cores to use to parallelize this function.
  #
  # Returns:  
  #   output: Either Primary, Secondary, or all FBP outputs in a data.frame
  #
  #############################################################################
  #hack to avoid Note about no visible binding for global variable ID
  #http://stackoverflow.com/questions/9439256/how-can-i-handle-r-cmd-check-no-visible-binding-for-global-variable-notes-when
  #look at the globalvariables() option or others in place of this issue
  # do not remove this comment until resolved
  ID <- NULL 
  if (!is.na(charmatch("input", search()))) {
    detach(input)
  }
  #If input is not provided, then calculate FBP with default values
  if (is.null(input)){
    fullList <- .FBPcalc(input)

  } else {
    #determine optimal number of pixels to process at each iteration
    if (is.null(m)){
      m <- ifelse(nrow(input) > 500000, 3000, 1000)
    }
    m <- ifelse(nrow(input) >= m, m, nrow(input))
    n0 <- round(nrow(input) / m)
    n <- ifelse(m * n0 >= nrow(input), n0, n0 + 1)
    #Set up parallel processing, if # of cores is entered
    if (cores > 1){
      #create and register a set of parallel R instances for foreach
      cl <- makeCluster(cores)
      registerDoParallel(cl)
      #process in parallel
      ca <- foreach(i=1:n, .packages='cffdrs') %dopar% {
        if (i==n){
          #Run FBP functions
          to.ls<-.FBPcalc(input[((i-1)*m+1):nrow(input),],output=output)
        }else {
          #Run FBP functions
          to.ls<-.FBPcalc(input[((i-1)*m+1):(i*m),],output=output)
        }
        to.ls
      }
      #close the processes
      stopCluster(cl)
      registerDoSEQ()
    #Run only a single process
    } else {
      ca <- vector('list',n)
      
      for(i in 1:n) {
        if(i == n) {
          foo <- input[((i - 1) * m + 1):nrow(input), ]
        } else {
          foo <- input[((i - 1) * m + 1):(i * m), ]
        }
        #Run FBP functions
        ca[[i]] <- .FBPcalc(foo, output = output)
      }    
    }
    #create a single keyed data table
    fullList <- rbindlist(ca)
    setkey(fullList, ID)
    #convert to data frame
    fullList <- as.data.frame(fullList)
  }
  
  return(fullList)
}

