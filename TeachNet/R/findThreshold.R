find.Threshold<- function(obs, stepsize=0.1, predict) {
  # finds best threshold for p
  
  # obs is vector of observations
  # stepsize for threshold
  # pred is vector of predictions
  
  # initialize values
  steps <- 1/stepsize + 1
  value <- rep(0, steps)
  threshold <- rep(0, steps)
  
  # calculate values
  for ( i in 1:steps) {
    threshold[i] = (i-1)*stepsize
    confusion <- confusion(predict,obs, threshold=threshold[i])
    value[i] = confusion[3,2]-confusion[3,1]
  }
  
  # structure outcome
  outcome <- structure(list(Value = value, Threshold = threshold), .Names = c("Value","Threshold"), row.names = c(NA,length(value)), class = "data.frame")
  
  #find maximum value
  maxi <- max(value)
  
  # select best threshold
  result <- subset(outcome, outcome$Value==maxi)
  
  return(result)
    
}
  