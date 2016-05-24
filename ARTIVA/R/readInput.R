readInput <-
function(data){
  # if data is character, i.e. a file path,
  if(is.character(data)){
    # if file exists
    if(file.exists(data)){
       DATA = as.matrix(read.table(data))
       return(DATA)
     } else {
       stop(paste("input file path",data,"is not valid\n"))
     }
  }else{
  # if data is already a matrix
  #if(is.matrix(data)){
    return(as.matrix(data))
  }
  # otherwise input data is incorrect :
  stop(paste("input data",data,"is not valid\n"))
}
