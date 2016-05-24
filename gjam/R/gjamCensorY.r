
gjamCensorY <- function(values,intervals,y,whichcol=c(1:ncol(y))){     
  
  #values    - in data set that are censored
  #intervals - matrix with 2 rows for lower and upper bounds for intervals
  
  censor <- vector('list',2)
  names(censor) <- c('columns','partition')
  
  censor$columns   <- whichcol
  censor$partition <- rbind(values,intervals)
  
  yc     <- y[,whichcol]
  
  for(k in 1:length(values)){
    yc[yc > intervals[1,k] & yc <= intervals[2,k]] <- values[k]  # censored data
  }
  
  list(censor = censor, y = yc)
}
