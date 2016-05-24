func.externalValidationAnalysis <-
function( input, output, result ){
  tmp <- NULL

  data <- data.frame(
    cbind(
      "obs_mean" = numeric(0),
      "obs" = numeric(0), 
      "pred" = numeric(0)
    )
  )

  tmp$colnames <- c( "observed mean", "observed value", "predicted value" )
  data.col <- matrix(                                                                                  
    c(paste( "C", 1:NCOL(data), sep =""), letters[1:NCOL(data)]),
    nrow = NCOL(data),
    ncol = 2,
    byrow = FALSE,
    dimnames = list(NULL, c("abrev", "trans"))
  )

  colnames(data) <- data.col[,"abrev"]
  data.col[1:NROW(tmp$colnames),"trans"] <- tmp$colnames

  tmp$pred$obs_col <- data.col[data.col[,"trans"]=="observed value", "abrev"]
  tmp$pred$pred_col <- data.col[data.col[,"trans"]=="predicted value", "abrev"]
  tmp$pred$obs_mean_col <- data.col[data.col[,"trans"]=="observed mean", "abrev"]
  
  tmp$pred$rows <- 1:NROW(input$predictData)
  tmp$pred$observed_mean <- mean(input$predictData[, ncol(input$predictData)])
  tmp$model$observed_mean <- mean(input$modelData[, ncol(input$modelData)])
  
  data[tmp$pred$rows, tmp$pred$pred_col] <- predict( result$fit$model, input$predictData )
  data[tmp$pred$rows, tmp$pred$obs_col] <- input$predictData[, ncol(input$predictData)]
  data[tmp$pred$rows, tmp$pred$obs_mean_col] <- rep(tmp$model$observed_mean, NROW(input$predictData))
  
  if( !is.null(output$writeTarget) ){
    writeLines("-- Start PARAMETER Prediction External Data Set --", con = output$writeTarget)
    writeLines("Data Table: ", con = output$writeTarget)
    cat( data.col[,"abrev"], "\n", sep="\t", file = output$writeTarget )
  
    write.table( round(data,output$round), file = output$writeTarget, sep="\t", row.names = FALSE, col.names = FALSE )
    writeLines("", con = output$writeTarget)
  
    writeLines("Data Table column names explanation:", con = output$writeTarget)
    write.table( data.col, file = output$writeTarget, sep="\t", row.names = FALSE, col.names=FALSE )
  
    writeLines("-- End PARAMETER Prediction External Data Set --", con = output$writeTarget)
    writeLines("", con = output$writeTarget)
  }
  
  tmp$data.stat <- func.get_data_stats( data, tmp$pred, "prediction", input$nu )
    
  #number of different test sets is missing
  return(
    list(
      "nTestSet" = input$nTestSet,
      "nTrainingSet" = input$nTrainingSet,
      
      #copy this values from stats
      "q2" = tmp$data.stat$x2, #q2 with y_mean(observed)
      "rmse" = tmp$data.stat$rmse,
      "observed_mean" = tmp$data.stat$observed_mean,
      "predicted_mean" = tmp$data.stat$predicted_mean,
      #rewrite this value(s)
      "nu" = input$nu,
      "data" = data,
      "data.col" = data.col,
      "TestSet" = input$predictData
    )
  )
}

