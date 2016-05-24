mainfunc.stats <-
function( data, obs = "observed", pred = "predicted", obs_mean = NULL, nu = 0, round = 4, stat_type, extOut = FALSE, extOutFile = NULL, func.call ){

  input <- list(
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
#  input <- func.completeInputPar( input, output, FALSE )
  result <- NULL
 
  if( !is.null(output$writeTarget) ){
    writeLines( "---- INPUT ----", con = output$writeTarget )
#    writeLines( "data ", con = output$writeTarget )
    write.table( data, file = output$writeTarget, sep="\t", row.names = FALSE )
    writeLines( "", con = output$writeTarget )
  }

  match.obs <- match(obs,colnames(data))
  match.pred <- match(pred,colnames(data))
  
  if( is.na(match.obs) && is.na(match.pred) ){
    cat("ERROR: no column",obs,"and no column",pred,"defined\n")
    cat("ERROR: assume that column #1 is",obs,"\n")
    cat("ERROR: assume that column #2 is",pred,"\n")
    colnames(data)[1]<-obs
    colnames(data)[2]<-pred
  }
  else{
    match.assume <- 1
    if( is.na(match.obs) ){
      cat("ERROR: no column",obs,"defined\n")
      if( match.pred == 1 )
        match.assume <- 2
      cat("ERROR: assume that column #",match.assume," is",obs,"\n")
      colnames(data)[match.assume]<-obs
    }
    if( is.na(match.pred) ){
      cat("ERROR: no column",pred,"defined\n")
      if( match.obs == 1 )
        match.assume <- 2
      cat("ERROR: assume that column #",match.assume," is",pred,"\n")
      colnames(data)[match.assume]<-pred
    }
  }  
  
  data.new <- data.frame(
    cbind(
      1:NROW(data), 
      data[,obs], 
      data[,pred]
    )
  )

  data.col <- matrix(
    c(paste( "C", 1:NCOL(data.new), sep =""), c("observed mean", "observed value", "predicted value") ),
    nrow = NCOL(data.new),
    ncol = 2,
    byrow = FALSE,
    dimnames = list(NULL, c("abrev", "trans"))
  )
  colnames(data.new) <- data.col[,"abrev"]
  
  tmp <- NULL
  tmp$obs_mean_col <- data.col[data.col[,"trans"]=="observed mean", "abrev"]
  tmp$obs_col <- data.col[data.col[,"trans"]=="observed value", "abrev"]
  tmp$pred_col <- data.col[data.col[,"trans"]=="predicted value", "abrev"]
  
  #set it in general
  if( is.null(obs_mean) ){ 
    #use the mean of all observed values
    data.new[,tmp$obs_mean_col] <- rep(mean(data.new[,tmp$obs_col]),NROW(data.new)) 
  }
  # not for stat_type == calibration
  if( stat_type == "prediction" ){
    if( is.character(obs_mean) ){ 
      #column is defined in data
      data.new[,tmp$obs_mean_col] = data[,obs_mean]
    }
    if( NROW(data) == NROW(obs_mean) ){ 
      #column is defined as list/vector
      data.new[,tmp$obs_mean_col] = obs_mean
    }
  }

  tmp$data.stat <- func.get_data_stats( data.new, tmp, stat_type, input$nu )
  
  if( stat_type == "calibration" ){
    result$fit <- 
      list(  
        #copy this values from stats
        "r2" = tmp$data.stat$x2,
        "rmse" = tmp$data.stat$rmse,
        "observed_mean" = tmp$data.stat$observed_mean,
        "predicted_mean" = tmp$data.stat$predicted_mean,
  
        "nu" = input$nu,
        "n" = nrow(data.new),
        "data.col" = data.col,
        "data" = data.new
      )      
  }

  if( stat_type == "prediction" ){
    result$pred <- 
      list(  
        #copy this values from stats
        "q2" = tmp$data.stat$x2,
        "rmse" = tmp$data.stat$rmse,
        "observed_mean" = tmp$data.stat$observed_mean,
        "predicted_mean" = tmp$data.stat$predicted_mean,
  
        "nu" = input$nu,
        "n" = nrow(data.new),
        "data.col" = data.col,
        "data" = data.new
      )      
  }

  if( output$toFile ){
    func.output.performanceValues( result, output )
    close(output$writeTarget)
  }
  
  # return as class and write it to stdout
  output$writeTarget = stdout()

  return( new("q2", result=result, output=output) )
}

