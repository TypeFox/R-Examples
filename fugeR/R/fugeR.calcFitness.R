# This function calculate the fitness of a fuzzy system
# According to the given weight for each criteria
#
# SENSITIVITY
# SPECIFICITY
# ACCURACY
# THRESHOLD
#
# RMSE
fugeR.calcFitness <-
function(predicted, actual, fuzzySystem) {
    #take fitness parameter and last best system found
    fugeR.lstFitParameter <- fugeRglobal[["fugeR.lstFitParameter"]]
    fugeR.BEST_FIT <- fugeRglobal[["fugeR.BEST_FIT"]]
    fugeR.BEST_SYS <- fugeRglobal[["fugeR.BEST_SYS"]]

    #get the verbose option
    verbose <- fugeRglobal[["fugeR.verbose"]]
    
    nbOut <- ncol(predicted)
    fitness <- c() #fitness to return
    rmse <- c() # Root mean square error
    accu <- c() # Accuracy
    sensi <- c() # Sensitivity
    speci <- c() # Specificity
    nbCase <- nrow(predicted)
    
    for(i in 1:nbOut) {
      #Classification
      AccuW  = fugeR.lstFitParameter$ACCU_W[i]
      sensiW = fugeR.lstFitParameter$SENSI_W[i]
      speciW = fugeR.lstFitParameter$SPECI_W[i]
      rmseW  = fugeR.lstFitParameter$RMSE_W[i]

      fitNum <- 0.0
      fitDenum <- 0.0
      #COMPUTE REGRESSION FIRST
      #Current out 
      currentPred <- predicted[,i]
      currentActual <- actual[,i]
      error <- (currentPred - currentActual)*(currentPred - currentActual)
      sumError <- sum(error)
      if( is.na(sumError) ) {
        stop("CHECK DATASET PLEASE, seems like all values are not numeric")
      }
      #COMPUTE RMSE
      mse <- sumError / nbCase
      rmse[i] <- sqrt(mse)
      fitNum <- fitNum + (rmse[i] * rmseW)
      fitDenum <- fitDenum + rmseW
      
      #COMPUTE CLASSIFICATION
      thresh <- fugeR.lstFitParameter$THRESH[i]
      
      if(!(is.na(thresh))) {
        perf <- cbinaryPerformance(currentPred, currentActual, thresh)
        
        accu[i]  <- perf$accu  # Accuracy
        sensi[i] <- perf$sensi # Sensitivity
        speci[i] <- perf$speci # Specificity
        
        fitNum <- fitNum + ((1-accu[i]) * AccuW)
        fitNum <- fitNum + ((1-sensi[i]) * sensiW)
        fitNum <- fitNum + ((1-speci[i]) * speciW)
        fitDenum <- fitDenum + AccuW
        fitDenum <- fitDenum + sensiW
        fitDenum <- fitDenum + speciW
      }
      fitness[i] <- fitNum / fitDenum
    }
    
    meanFitness <- mean(fitness)

    #KEEP THE BEST FUZZY
    if (meanFitness < fugeR.BEST_FIT | is.na(fugeR.BEST_FIT)) {
      fugeR.BEST_SYS <- fuzzySystem
      fugeR.BEST_FIT <- meanFitness
      fugeR.BEST_SYS$fitness <- fitness
      fugeR.BEST_SYS$rmse <- rmse
      fugeR.BEST_SYS$accu <- accu
      fugeR.BEST_SYS$sensi <- sensi
      fugeR.BEST_SYS$speci <- speci
      
      assign('fugeR.BEST_SYS', fugeR.BEST_SYS, envir=fugeRglobal)
      assign('fugeR.BEST_FIT', fugeR.BEST_FIT, envir=fugeRglobal)
    
      if(verbose) {
       cat('RMSE : \t',rmse, '\n')
        if(!(is.na(thresh))) {
          cat('Sensitivity : \t',sensi, '\n')
          cat('Specificity : \t',speci, '\n')
        }
      }
    }
    
   #RETURN the mean fitness of each output variable
   return(meanFitness)
} 
