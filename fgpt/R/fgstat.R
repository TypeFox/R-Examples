fgstat <-
function(rand,marks,FUN=mean, ...){
  if(missing(marks)){
    output <- sapply(rand, FUN, ...)  
  }else if(is.vector(marks)){
    output <- sapply(rand, function(x){FUN(marks[x], ...)})
  }else if(dim(marks)[2]==2){
    output <- sapply(rand, function(x){FUN(marks[,1],marks[x,2], ...)}) 
  }else{
    output <- sapply(rand, function(x){FUN(diag(marks[x,1:dim(marks)[1]]), ...)}) 
  }
  if(is.na(sum(output))){
    warning("This function produced NAs." )
  }
  return(output)
}
