func.crossValidationAnalysis <-
function( input, output ){
  tmp <- NULL
  samp <- NULL
  
  cv.data <- data.frame(
    cbind(
      "no" = numeric(0),
      "obs_mean" = numeric(0), 
      "obs" = numeric(0), 
      "pred" = numeric(0)
    )
  )

  # extend the dataframe for the regression values
  cv.data <- func.extendDataframeForRegValues( cv.data, output$coefficients )
  cv.data.col <- matrix(
    c(paste( "C", 1:NCOL(cv.data), sep =""), letters[1:NCOL(cv.data)]),
    nrow = NCOL(cv.data),
    ncol = 2,
    byrow = FALSE,
    dimnames = list(NULL, c("abrev", "trans"))
  )
  colnames(cv.data) <- cv.data.col[,"abrev"]

  tmp$colnames <- c( "no in modelData", "observed mean", "observed value", "predicted value" )

  cv.data.col[1:NROW(tmp$colnames),"trans"] <- tmp$colnames
  tmp$colnames.i <- NROW(tmp$colnames) + 1

  if(names(output$coefficients)[1] == "(Intercept)"){
    cv.data.col[tmp$colnames.i,"trans"] <- "const"
    increment(tmp$colnames.i)
  }
  
  for( i in tmp$colnames.i:NCOL(cv.data) )
    cv.data.col[i,"trans"] <- letters[i - tmp$colnames.i + 1]

  tmp$cv$no_col <- cv.data.col[cv.data.col[,"trans"]=="no in modelData", "abrev"]
  tmp$cv$pred_col <- cv.data.col[cv.data.col[,"trans"]=="predicted value", "abrev"]
  tmp$cv$obs_col <- cv.data.col[cv.data.col[,"trans"]=="observed value", "abrev"]
  tmp$cv$obs_mean_col <- cv.data.col[cv.data.col[,"trans"]=="observed mean", "abrev"]

  #perform cross validation
  samp$no = 1 
  input$sampleOrder <- NULL

  for( iRun in 1:input$nRun ){
    samp$orderThisRun <- func.createSample( input ) 
    input$runSampleOrder[iRun] <- list( samp$orderThisRun )
  
    while(NROW(samp$orderThisRun) > 0){
      if( NROW(samp$orderThisRun)%%input$nTestSet == 0)
        samp$rows <- samp$orderThisRun[c(1:input$nTestSet)]
      else
        samp$rows <- samp$orderThisRun[c(1:(input$nTestSet - 1))]
    
      samp$trainingSet <- input$modelData[-samp$rows,] #regressionSet
      samp$testSet <- input$modelData[samp$rows[order(samp$rows)],]  #predictionSet, ordered

      tmp$testSet[samp$no] <- list( samp$testSet )
  
      samp$model <- glm( formula = input$regressionFormula, data=samp$trainingSet )
      
      if(input$nRun == 1)
        tmp$rowNo <- rownames(samp$testSet)
      else
        tmp$rowNo <- 1:nrow(samp$testSet)+NROW(cv.data)

      cv.data[tmp$rowNo, tmp$cv$pred_col] <- predict(samp$model, samp$testSet)
      cv.data[tmp$rowNo, tmp$cv$obs_col] <- samp$testSet[, ncol(input$modelData)]
      cv.data[tmp$rowNo, tmp$cv$obs_mean_col] <- rep( mean(samp$trainingSet[, ncol(input$modelData)]), NROW(tmp$rowNo))
      cv.data[tmp$rowNo, tmp$cv$no_col] <- as.numeric(rownames(samp$testSet))
     
      for(c in 1:NROW(output$coefficients) )
        cv.data[tmp$rowNo, ncol(cv.data)-NROW(output$coefficients)+c] <- rep( samp$model$coefficients[c], NROW(tmp$rowNo))
                                                                                               
      samp$orderThisRun <- samp$orderThisRun[-c(1:NROW(samp$rows))]
  
      if( !is.null(output$writeTarget) ){
        cat("Sample",samp$no,":",sort(samp$rows),"\n", file = output$writeTarget)
        cat("Leave-X-Out Mean",mean(samp$trainingSet[,ncol(input$modelData)]),"\n", file = output$writeTarget)
        func.output.regressionFormulaWithCoefficients(output, colnames(input$modelData))

        writeLines("", con = output$writeTarget)
      }
      increment(samp$no)
    }

    increment(iRun)
  }

  #copy it unsorted  
  pred.data <- cv.data
#  print(order(pred.data[tmp$cv$no_col,]))
  pred.data <- pred.data[order(pred.data[,tmp$cv$no_col]),]
  pred.data <- pred.data[ ,c(tmp$cv$no_col, tmp$cv$obs_col, tmp$cv$pred_col) ]

  # sort the cv.data in case of sample fractioning
#  pred.data <- pred.data[order(pred.data[tmp$cv$no_col,]),]
#  print(cv.data)
#  print(pred.data)

  tmp$pred$colnames <- c( "no in modelData", "observed value", "predicted value" )

  pred.data.col <- matrix(                                                                                  
    c(paste( "C", 1:NCOL(pred.data), sep =""), letters[1:NCOL(pred.data)]),
    nrow = NCOL(pred.data),
    ncol = 2,
    byrow = FALSE,
    dimnames = list(NULL, c("abrev", "trans"))
  )

  colnames(pred.data) <- pred.data.col[,"abrev"]
  pred.data.col[1:NROW(tmp$pred$colnames),"trans"] <- tmp$pred$colnames

  if( !is.null(output$writeTarget) ){
    writeLines("-- Start OVERVIEW Parameter Cross Validation --", con = output$writeTarget)
    writeLines("", con = output$writeTarget)
    #use coefficients from glm, should be same for cross valdiation
    func.output.linearFormula( output )
    writeLines("", con = output$writeTarget)
  }
#  print(output$coefficients)

  if( !is.null(output$writeTarget) ){
    writeLines("Data Table: ", con = output$writeTarget)
    cat( cv.data.col[,"abrev"], "\n", sep="\t", file = output$writeTarget )
  
    write.table( round(cv.data,output$round), file = output$writeTarget, sep="\t", row.names = FALSE, col.names = FALSE )
    writeLines("", con = output$writeTarget)
  
    writeLines("Data Table column names explanation:", con = output$writeTarget)
    write.table( cv.data.col, file = output$writeTarget, sep="\t", row.names = FALSE, col.names=FALSE )
  
    writeLines("-- End OVERVIEW Parameter Cross Validation --", con = output$writeTarget)
    writeLines("", con = output$writeTarget)
  }
  
  #hier koennten theoretisch noch die einzelnen samples und ihre ergebnisse rein
  tmp$res$cv = list(
    #redundant to res$pred
    "nTestSet" = input$nTestSet,
    "nTrainingSet" = input$nTrainingSet,
    #
    "nFold" = input$nFold,
    "testSetSizeVaries" = input$testSetSizeVaries,
    #rewrite this value(s)
    "nRun" = input$nRun,
    
    "data" = cv.data,
    "data.col" = cv.data.col,
    "TestSet" = tmp$testSet
  )
  
  tmp$data.stat <- func.get_data_stats( cv.data, tmp$cv, "prediction", input$nu )
    
  tmp$res$pred = list(
    "nTestSet" = input$nTestSet,
    "nTrainingSet" = input$nTrainingSet,
    #copy this values from stats
    "q2" = tmp$data.stat$x2, #q2 with y_mean_i-1(training) - cv 
    "rmse" = tmp$data.stat$rmse,
    "observed_mean" = tmp$data.stat$observed_mean,
    "predicted_mean" = tmp$data.stat$predicted_mean,
    "nu" = input$nu,
    "data" = pred.data,
    "data.col" = pred.data.col
  )

  #number of different test sets is missing
  return( tmp$res )
}

