insertRow <-
function(dataframe, newrow, index){
  if(!(index == nrow(dataframe)+1)){
    dataframe[seq(index+1,(nrow(dataframe)+1)),] <- dataframe[seq(index,nrow(dataframe)),]
  }
  dataframe[index,] <- newrow
  return(dataframe)
}
