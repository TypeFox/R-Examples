func.output.performanceValues <-
function(result, output){
  varyTrainingSet = ""
  varyTestSet = ""

  writeLines("\n---- CALL ----", con = output$writeTarget)
  writeLines(deparse(output$call), con = output$writeTarget)

  writeLines("\n---- RESULTS ----\n", con = output$writeTarget)

  if( !is.null(result$fit) ){
    writeLines( "-- MODEL CALIBRATION (linear regression)", con = output$writeTarget )
    cat( "#Elements: \t",result$fit$n,"\n\n",sep="", file = output$writeTarget )
    
    nFormat <- format(c("observed_mean"=result$fit$observed_mean, "predicted_mean"=result$fit$predicted_mean, "rmse"=result$fit$rmse, "r2"=result$fit$r2), digit=output$round, nsmall = output$round)
    
    cat("mean (observed): \t", nFormat["observed_mean"], "\n", sep="", file = output$writeTarget )
    cat("mean (predicted): \t", nFormat["predicted_mean"], "\n", sep="", file = output$writeTarget )
    cat("rmse (nu = ",result$fit$nu,"): \t\t", nFormat["rmse"], "\n", sep="", file = output$writeTarget )
    cat("r^2: \t\t\t", nFormat["r2"], "\n", sep="", file = output$writeTarget )
  }
  
  out.perf = "model and prediction set available"
  if( !is.null(result$cv) )
    out.perf = "cross validation"
  
  if( !is.null(result$pred) || !is.null(result$cv) )
    cat( "\n-- PREDICTION PERFORMANCE (", out.perf, ")\n", sep="", file = output$writeTarget )
  
  if( !is.null(result$cv) ){
    #cross validation specific values
    cat("#Runs: \t\t\t\t", result$cv$nRun, "\n", sep="", file = output$writeTarget )
    cat("#Groups: \t\t\t", result$cv$nFold, "\n", sep="", file = output$writeTarget )
    
#    y_means der gruppen
  
    if( result$cv$testSetSizeVaries ){
      varyTrainingSet = " (+1)"
      varyTestSet = " (-1)"
    }
  }
  
  out.train = "#Elements Model Set: \t\t"
  out.test = "#Elements Prediction Set: \t" 
  
  if( !is.null(result$cv) ){
    out.train = "#Elements Training Set: \t"
    out.test = "#Elements Test Set: \t\t"
  }

  if( !is.null(result$pred) ){
    if( !is.null(result$pred$nTrainingSet) && !is.null(result$pred$nTrainingSet) ){
      cat( out.train, result$pred$nTrainingSet,varyTrainingSet, "\n",sep="", file = output$writeTarget )
      cat( out.test, result$pred$nTestSet,varyTestSet, "\n\n",sep="", file = output$writeTarget )
    }
  
    nFormat <- format(c("observed_mean"=result$pred$observed_mean, "predicted_mean"=result$pred$predicted_mean, "rmse"=result$pred$rmse, "q2"=result$pred$q2), digits = output$round, nsmall = output$round)
    
    cat("mean (observed): \t", nFormat["observed_mean"], "\n", sep="", file = output$writeTarget )
    cat("mean (predicted): \t", nFormat["predicted_mean"], "\n", sep="", file = output$writeTarget )
    cat("rmse (nu = ",result$pred$nu,"): \t\t", nFormat["rmse"], "\n",sep="", file = output$writeTarget )

    # calculated with Y_mean^training
    cat("q^2: \t\t\t", nFormat["q2"], "\n",sep="", file = output$writeTarget )
  }
}

