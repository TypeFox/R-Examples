#ClassificationResults class
  #Calculate and stores the results of a regression algorithm

RegressionResults <- R6::R6Class("RegressionResults",

  public = list(

    #Real and predicted classes
    predictions = NULL,

    #Metrics
    MAE = NaN,
    MSE = NaN,
    RMSE = NaN,
    R2 = NaN,
    #add more metrics...

    initialize = function(predictions){

      self$predictions <- predictions

      #Calculate metrics
      self$calculateMetrics()
    },


    #Calculate metrics
    calculateMetrics = function(){

      MAE <- 0
      MSE <- 0
      RMSE <- 0
      R2 <- 0

      actualMean <- mean(as.numeric(as.character(self$predictions$Real)))
      predictedMean <- mean(as.numeric(as.character(self$predictions$Predicted)))

      SSR <- 0
      SST <- 0

      actual_i <- NaN
      predicted_i <- NaN

      for(i in 1:nrow(self$predictions)){
        actual_i <- as.numeric(as.character(self$predictions[i,1]))
        predicted_i <- as.numeric(as.character(self$predictions[i,2]))

        MAE <- MAE + abs(predicted_i - actual_i) #predicted - actual
        MSE <- MSE + (predicted_i - actual_i)**2
        SST <- SST + (actual_i - actualMean)**2
        SSR <- SSR + (predicted_i - predictedMean)**2
      }

      MAE <- MAE/nrow(self$predictions)
      self$MAE <- MAE
      MSE <- MSE/nrow(self$predictions)
      self$MSE <- MSE
      self$RMSE <- sqrt(MSE)
      self$R2 <- SSR/SST
    },

    #Print object
    print = function() {

      cat("MAE", ": ", self$MAE, "\n", sep="")
      cat("MSE", ": ", self$MSE, "\n", sep="")
      cat("RMSE", ": ", self$RMSE, "\n", sep="")
      cat("R2", ": ", self$R2, "\n", sep="")
    },

    #Export results
    export = function(filename) {

      text <- ""
      text <- paste0(text, "MAE", ": ", self$MAE, "\n")
      text <- paste0(text, "MSE", ": ", self$MSE, "\n")
      text <- paste0(text, "RMSE", ": ", self$RMSE, "\n")
      text <- paste0(text, "R2", ": ", self$R2, "\n")
      text <- paste0(text, "\n")

      #Save results
      fileConn<-file(filename)
      writeLines(text, fileConn)
      close(fileConn)

    }

  ),

  private = list(

  )
)
