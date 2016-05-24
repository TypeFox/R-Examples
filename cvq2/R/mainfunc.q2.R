mainfunc.q2 <-
function( modelData, predictData = NULL, formula = NULL, nFold = N, nRun = 1,
nu = 0, round = 4, extOut = FALSE, extOutFile = NULL, func.call ){
#  call <- match.call()

  apply_cv <- FALSE
  if( is.null(predictData) ){
    apply_cv <- TRUE
    predictData <- modelData
  }

  N <- nrow(predictData)
  
  input <- list(
    "modelData" = modelData,
    "predictData" = predictData,
    "regressionFormula" = formula,
    "nFold" = nFold,
    "nRun" = nRun,
    #degrees of freedom
    "nu" = nu
  )

  output <- list(
    "call" = func.call,
    "round" = round,
    "toFile" = FALSE,
    "writeTarget" = NULL,
    "coefficients" = NULL
  )

  output <- func.completeOutputPar( output, extOut, extOutFile )
  input <- func.completeInputPar( input, output, apply_cv )
 
  result <- NULL
  tmp <- NULL

  if( !is.null(output$writeTarget) ){
    writeLines( "---- INPUT ----", con = output$writeTarget )
    writeLines( "Model Data: ", con = output$writeTarget )
    write.table( input$modelData, file = output$writeTarget, sep="\t", row.names = FALSE )

    if( !apply_cv ){
      writeLines( "Prediction Data: ", con = output$writeTarget )
      write.table( input$predictData, file = output$writeTarget, sep="\t", row.names = FALSE )
    }
    writeLines( "", con = output$writeTarget )
  }

  # linear regression
  result$fit <- func.linearRegressionAnalysis( input, output )

  #the cross validation, as predictDataSet is missing
  if( apply_cv ){
    input$predictData <- NULL

    # leave-X-out, cross validation
    #receive $cv and $pred
    tmp.result <- func.crossValidationAnalysis( input, output )
    result$cv <- tmp.result$cv
    result$pred <- tmp.result$pred
  }
  # the actual prediction
  else{
    # use fit from linear regression, validate it with external data set
    result$pred <- func.externalValidationAnalysis( input, output, result )
  }

  if( output$toFile ){
    func.output.performanceValues( result, output )
    close(output$writeTarget)
  }
  
  # return as class and write it to stdout
  output$writeTarget = stdout()

  if( apply_cv )
    return( new("cvq2", result=result, output=output) )
  else
    return( new("q2", result=result, output=output) )
}

