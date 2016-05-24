#ClassificationResults class
  #Calculate and stores the results of a classification algorithm

ClassificationResults <- R6::R6Class("ClassificationResults",

  public = list(

    #Real and predicted classes
    predictions = NULL,

    #confusion matrix
    confusionMatrix = NULL,

    #Num of unclassified instances
    unclassified = 0,

    #Metrics
    accuracy = NaN,
    precision = NaN,
    recall = NaN,
    FMeasure = NaN,
    #add more metrics...

    initialize = function(predictions){

      self$predictions <- predictions

      #Compute confusion matrix
      self$calculateConfusionMatrix()

      #Calculate metrics
      self$calculateMetrics()
    },


    #Calculate confusion matrix
    calculateConfusionMatrix = function(){

      #Join real and predicted classes to ensure that are all
      classes <- unique(c(self$predictions$Real, self$predictions$Predicted))
      numClasses <- length(classes)

      #
      #   | PREDICTED CLASSES
      # __|___________________
      #   |
      # R |
      # E |    m a t r i x
      # A |    v a l u e s
      # L |
      #   |
      #

      #Create vector of zeros for each cell
      vector <- c()
      for(i in 1:(numClasses*numClasses))
        vector <- c(vector, 0)

      #Create matrix
      self$confusionMatrix <- matrix(vector, nrow=numClasses, ncol=numClasses)
      #Assign class names
      rownames(self$confusionMatrix) <- classes
      colnames(self$confusionMatrix) <- classes

      #Fill confusion matrix
      for(i in 1:length(self$predictions$Real)){

        if(nchar(self$predictions$Predicted[i]) > 0){
          if(tolower(self$predictions$Predicted[i]) == "unclassified") {
            self$unclassified <- self$unclassified + 1
          }
          else {
            realIndex <- -1
            predictedIndex <- -1
            for(c in 1:numClasses){
              if(self$predictions$Real[i] == classes[c]){
                realIndex <- c
              }
              if(self$predictions$Predicted[i] == classes[c]){
                predictedIndex <- c
              }
            }
            self$confusionMatrix[realIndex, predictedIndex] <- self$confusionMatrix[realIndex, predictedIndex]+1
          }
        }
      }

    },

    calculateMetrics = function(){

      #Confusion matrix must be calculated
      if(is.null(self$confusionMatrix)){
        self$calculateConfusionMatrix()
      }

      #Calculate different metrics
      self$accuracy <- private$calculateAccuracy()
      self$precision <- private$calculatePrecision()
      self$recall <- private$calculateRecall()
      self$FMeasure <- private$calculateFMeasure()

    },

    #Print object
    print = function() {

      cat("Confusion matrix: ", "\n", sep="")
      print(self$confusionMatrix)
      cat("\n", sep="")
      cat(" Unclassified istances: " , self$unclassified, "\n", sep="")
      cat("\n", sep="")

      cat("Metrics: ", "\n", sep="")
      cat("   ", "Accuracy", ": ", self$accuracy, "\n", sep="")
      cat("   ", "Precision", ": ", self$precision, "\n", sep="")
      cat("   ", "Recall", ": ", self$recall, "\n", sep="")
      cat("   ", "FMeasure", ": ", self$FMeasure, "\n", sep="")
    },

    #Export results
    export = function(filename) {

      text <- ""
      text <- paste0(text, "Metrics: ", "\n")
      text <- paste0(text, "   ", "Accuracy", ": ", self$accuracy, "\n")
      text <- paste0(text, "   ", "Precision", ": ", self$precision, "\n")
      text <- paste0(text, "   ", "Recall", ": ", self$recall, "\n")
      text <- paste0(text, "   ", "FMeasure", ": ", self$FMeasure, "\n")
      text <- paste0(text, "\n")

      #Save results
      fileConn<-file(filename)
      writeLines(text, fileConn)
      close(fileConn)

    }

  ),

  private = list(

    #Calculate accuracy
    calculateAccuracy = function() {

      if(sum(self$confusionMatrix) > 0) {
        return(sum(diag(self$confusionMatrix))/sum(self$confusionMatrix))
      }
      else {
        return(NaN)
      }
    },

    #Calculate precision
    calculatePrecision = function() {

      #If there are two classes, calculate for positive class
      if(ncol(self$confusionMatrix) == 2){
        precision <- self$confusionMatrix[1, 1] / sum(self$confusionMatrix[,1])
      }
      #If there are more than two, calculate for each class
      else{
        precision <- 0
        n <- 0

        for(i in 1:ncol(self$confusionMatrix)){
          if(sum(self$confusionMatrix[,i]) != 0){
              precision <- precision + (self$confusionMatrix[i, i] / sum(self$confusionMatrix[,i]))
              n <- n+1
          }
        }

        if(n > 0){
          precision <- precision / n
        }
        else{
          precision <- NaN
        }
      }

      return(precision)
    },

    #Calculate recall
    calculateRecall = function() {

      #If there are two classes, calculate for positive class
      if(ncol(self$confusionMatrix) == 2){
        recall <- self$confusionMatrix[1, 1] / sum(self$confusionMatrix[1,])
      }
      #If there are more than two, calculate for each class
      else{
        recall <- 0
        n <- 0

        for(i in 1:ncol(self$confusionMatrix)){
          if(sum(self$confusionMatrix[i,]) != 0){
            recall <- recall + (self$confusionMatrix[i, i] / sum(self$confusionMatrix[i,]))
            n <- n+1
          }
        }

        if(n > 0) {
          recall <- recall/n
        }
        else {
          recall <- NaN
        }
      }

      return(recall)
    },

    #Calculate FMeasure
    calculateFMeasure = function() {

      FMeasure <- 0

      precision <- private$calculatePrecision()
      recall <- private$calculateRecall()

      if((precision + recall) > 0) {
        FMeasure <- 2 * ((precision*recall)/(precision+recall))
      }
      else {
        FMeasure <- NaN
      }

      return(FMeasure)
    }

  )
)
