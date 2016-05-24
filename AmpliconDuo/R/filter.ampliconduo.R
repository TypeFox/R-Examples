filter.ampliconduo <-
function(x, min.freq = 1, OR = NULL, q = NULL, p = NULL, remove = FALSE){
  data <- x
  # x is a data.frame with the ampliconduo data
  filter <- which(data$freqA < min.freq | data$freqB < min.freq)
  data$rejected[filter]<- TRUE
  
  if(!is.null(OR)){
    filter <- which(data$OR < OR)
    data$rejected[filter]<- TRUE
  }
  
  if(!is.null(q)){
    filter <- which(data$q < q)
    data$rejected[filter]<- TRUE
  } 
  
  if(!is.null(q)){
    filter <- which(data$p < p)
    data$rejected[filter]<- TRUE
  } 
  
  if(remove){
    data <- data[data$rejected == FALSE, ]
  }
  
  return(data)
}
