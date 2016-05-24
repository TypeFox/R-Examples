estimationFunction <-
function(daten, decs = 4) {
  if(dim(daten)[1]!=dim(daten)[2]) stop("Number of rows does not equal number of columns")
  datenMatrix = as.matrix(daten)[, 1:length(daten)]
  datenMatrix = matrix(datenMatrix, length(daten), length(daten), byrow = T)
  
  sMatrix = createSmallerMatrix(datenMatrix)
  sMatrix = fillMatrix(datenMatrix, sMatrix)
  
  nrOfZerosSmaller = (length(which(sMatrix == 0)))*9 ## Number of zeros overall in bigger matrix
  nrOfZerosBigger = as.numeric(length(which(datenMatrix == 0)))
  nrBoxes = (length(which(sMatrix == 1)))*9
  
  ## Calculating the estimator for q
  qEst = max(round((nrOfZerosBigger - nrOfZerosSmaller)/nrBoxes, decs), 0)
  
  sMatrix2 = createSmallerMatrix(sMatrix)
  sMatrix2 = fillMatrix(sMatrix, sMatrix2)
  
  nrOfZerosSmaller = (length(which(sMatrix2 == 0)))*9 ## Number of zeros overall in bigger matrix
  nrOfZerosBigger = as.numeric(length(which(sMatrix == 0)))
  nrBoxes = (length(which(sMatrix2 == 1)))*9
  
  ## Calculating the estimator for p
  pEst = max(round((nrOfZerosBigger - nrOfZerosSmaller)/nrBoxes, decs), 0)
  result = list(pEst, qEst)
  names(result) = c("p-Estimator", "q-Estimator")
  return(result)
}
