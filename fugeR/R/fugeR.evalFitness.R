# This function calculate the fitness of a fuzzy system
# According to the given weight for each criteria
#
# SENSITIVITY
# SPECIFICITY
# ACCURACY
# THRESHOLD
#
# RMSE
fugeR.evalFitness <-
    function(predicted, actual, fuzzySystem, lstParameter) {
        #take fitness parameter and last best system found
        #fugeR.lstFitParameter <- fugeRglobal[["fugeR.lstFitParameter"]]
        
        nbOut <- ncol(predicted)
        fitness <- c() #fitness to return
        rmse <- c() # Root mean square error
        accu <- c() # Accuracy
        sensi <- c() # Sensitivity
        speci <- c() # Specificity
        nbCase <- nrow(predicted)
        
        for(i in 1:nbOut) {
            #Classification
            AccuW  = lstParameter$ACCU_W[i]
            sensiW = lstParameter$SENSI_W[i]
            speciW = lstParameter$SPECI_W[i]
            rmseW  = lstParameter$RMSE_W[i]
            
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
            thresh <- lstParameter$THRESH[i]
            
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

        res <- list( meanFitness=meanFitness,
                     fitness=fitness,
                     rmse=rmse,
                     accu=accu,
                     sensi=sensi,
                     speci=speci
                   )
        
        #RETURN the mean fitness of each output variable
        return(res)
    } 
