sample2pop <- function(ctmatrix, population){
  
  if(is.null(colnames(ctmatrix)) | is.null(rownames(ctmatrix)))
    stop("colnames and rownames of ctmatrix must be defined as integer identifiers to match with first column of population")
  if ((dim(population)[2] != 2)) 
    stop("population must have two columns")
  if ((dim(ctmatrix)[1] != dim(ctmatrix)[2])) 
    stop("ctmatrix must be a square matrix")
  if ((dim(population)[1] != dim(ctmatrix)[1])) 
    stop("number of rows of population and ctmatrix must be equal")
  
  for (i in 1:dim(population)[1]){
    j <- which(colnames(ctmatrix) == population[i, 1])
    ctmatrix[j,] <- population[j, 2]*ctmatrix[j,]/sum(ctmatrix[j,])
  }
return(ctmatrix)  
}
