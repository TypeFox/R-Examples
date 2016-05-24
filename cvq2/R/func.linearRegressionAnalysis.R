func.linearRegressionAnalysis <-
function( input, output ){
  tmp <- NULL

  #data.meanY
  observed_mean = mean(input$modelData[,ncol(input$modelData)])

  model <- glm( formula = input$regressionFormula, data = input$modelData )
  #replace the coefficients outside this function too
  eval.parent(substitute(output$coefficients <- model$coefficients))

  if( !is.null(output$writeTarget) ){
    writeLines("\n---- Start linear regression ----", con = output$writeTarget)
    func.output.regressionFormulaWithCoefficients( output, colnames(input$modelData) ) 
    writeLines( "", con = output$writeTarget )
  }
  
  data.col <- matrix(
    c(paste( "C", 1:2, sep =""), c("observed value", "predicted value") ),
    nrow = 2,
    ncol = 2,
    byrow = FALSE,
    dimnames = list(NULL, c("abrev", "trans"))
  )
  
  data <- cbind( input$modelData[,NCOL(input$modelData)], 1:nrow(input$modelData) )
  colnames(data) <- data.col[,"abrev"]

  tmp$obs_col <- data.col[data.col[,"trans"]=="observed value", "abrev"]
  tmp$pred_col <- data.col[data.col[,"trans"]=="predicted value", "abrev"]
  
  data[, tmp$pred_col] <- predict( model, input$modelData )

  if( !is.null(output$writeTarget) ){
    writeLines("Observed vs. predicted values: ", con = output$writeTarget)
    cat( data.col[,"abrev"], "\n", sep="\t", file = output$writeTarget )
    write.table( round(data,output$round), file = output$writeTarget, sep="\t", row.names = FALSE, col.names = FALSE )
  
    writeLines("", con = output$writeTarget)
    writeLines("Table column names explanation:", con = output$writeTarget)
    write.table( data.col, file = output$writeTarget, sep="\t", row.names = FALSE, col.names=FALSE )

    writeLines("---- End linear regression ----", con = output$writeTarget)
    writeLines("", con = output$writeTarget)
  }
  
  #linear regression has per default a nu of 0
  tmp$nu = 0
  tmp$data.stat <- func.get_data_stats( data, tmp, "calibration", tmp$nu )
      
  return( 
    list(  
      #copy this values from stats
      "r2" = tmp$data.stat$x2,
      "rmse" = tmp$data.stat$rmse,
      "observed_mean" = tmp$data.stat$observed_mean,
      "predicted_mean" = tmp$data.stat$predicted_mean,

      "n" = nrow(input$modelData),
      "nu" = tmp$nu,
      "data.col" = data.col,
      "data" = data,
      "model" = model
    )      
  )
}

