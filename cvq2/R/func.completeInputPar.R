func.completeInputPar <-
function( input, output, apply_cv ){
#    input$splitSizeDataSet = floor( nrow(data) / input$nFold )
  input$testSetSizeVaries = FALSE 
 
  if( apply_cv ){
    # minimum 2 groups
    # maximum N groups == Leave one out cross validation
    # there can not be more groups than elements in the data set
    if( input$nFold < 2 || input$nFold > nrow(input$predictData) ){
      cat("It is not possible, to have less than 2 groups or more groups (",input$nFold,") than elements in the data set exist (",nrow(input$predictData),").\n")
      cat("In this case, 2 <= nFold <=",nrow(input$predictData),"\n")
      stop("Change parameter settings and start againg")
    }
    # maximum N groups == Leave one out cross validation
    # in every iteration (run) one get the same distribution of training and test set
    # therefore it is not necessary to perform more than one run
    if( input$nFold == nrow(input$predictData) && input$nRun > 1 ){
      cat("Training and test set distribution will be equal for every individual run. Therefore it is not necessary to perform more than one run, nRun is set to 1\n")
      input$nRun <- 1
    }
    
    if( input$nRun < 1 ){
      cat("nRun (",input$nRun,") must not be smaller than 1, set it to 1\n")
      input$nRun <- 1
    }
    
    input$nTestSet = ceiling( nrow(input$predictData) / input$nFold )
    input$nTrainingSet = nrow(input$modelData) - input$nTestSet
    # everytime the same size for training and test set
  
    nTrainingSetMin = input$nTrainingSet
  
    # no equal distribution size for test and training set possible
    # distribution can vary (test set + 1), (training set - 1)
    if( input$nTestSet != NROW(input$predictData) / input$nFold ){
      input$testSetSizeVaries = TRUE
      decrement(nTrainingSetMin)
    }
  }
  else{
    input$nTestSet = nrow(input$predictData)
    input$nTrainingSet = nrow(input$modelData)
 
    if( input$nRun != 1 ){
      cat("Use different prediction data compared to model data. Therefore it is not necessary to perform more than one run, nRun is set to 1\n")
      input$nRun <- 1
    }
  }

  # number of x + 1 is minimum for training set
  # here ncol == (number of x , y)
  if( input$nTrainingSet < ncol(input$modelData) ){
    cat("min(nTrainingSet) (",nTrainingSetMin,") is to small to create a linear model, must be at least (",ncol(input$modelData),")\n")
    stop("Change parameter settings and start againg")
  }

  #identify generic formula from data
  if( is.null(input$regressionFormula) )
    input$regressionFormula = func.constructRegressionFormula( colnames(input$modelData) )

  # exclude all lines with error
  input$modelData <- input$modelData[rowSums(is.na(input$modelData)) == 0,]

  # reorder the columns to match the given order in formula
  input$modelData <- func.sortDataColumns( input, input$modelData, output$writeTarget )
  input$predictData <- func.sortDataColumns( input, input$predictData, output$writeTarget )

  return( input )
}

