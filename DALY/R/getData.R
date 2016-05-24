## Retrieve user input from table
## Return input as matrix

getData <-
function(x, type){
  if (type == "pop") rArray <- array(0, dim = c(5, 2))
  if (type == "data") rArray <- array(0, dim = c(5, 6))

  if (!all(is.na(x))){
    for (i in seq(dim(x)[1]))
	  for (j in seq(dim(x)[2]))
	    rArray[i, j] <- ifelse(is.na(x[i, j]), 0, x[i, j])
  }

  return(rArray)
}