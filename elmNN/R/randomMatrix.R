randomMatrix <-
function(nCols,nRows,minValue,maxValue) {
  myMat<-matrix(runif(nCols*nRows,min=minValue,max=maxValue), ncol=nCols)
  myMat 
}
