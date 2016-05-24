getDataLongi <- function(theData, numTime, sizeSubsetGetData){
  #a list of data frames
  dataList <- list()


  #get the index of instances to draw
  the_index <- sample(1:(nrow(theData) / (numTime - 1)),
                      sizeSubsetGetData, replace=FALSE)

  #initialize data frames
  for (i in 1:(numTime - 1)){
    dataList[[i]] <- data.frame()
  }

  for (i in 1:length(the_index)){
    for (j in 1:length(dataList)){
      dataList[[j]] <- rbind(dataList[[j]],
                             theData[the_index[i] +
                                       ((j - 1) * (nrow(theData) /
                                                     (numTime - 1))), ])
    }
  }

  resData <- do.call(rbind, dataList)
  return(resData)
}


getDataCross <- function(theData, sizeSubset) {

  newData <- as.data.frame(theData[sample(1:nrow(theData),
                                          sizeSubset,
                                          replace=FALSE),])


  #check if any column has the same constant
  while (any(apply(newData, 2, sameValue))) {
    newData <- as.data.frame(theData[sample(1:nrow(theData),
                                            sizeSubset,
                                            replace=FALSE),])
  }

  return(newData)
}

sameValue <- function(x) {
  anySame <- FALSE

  #if all datas instance is the same
  if (all(x[1] == x[2:length(x)])) {
    anySame <- TRUE
  }

  return(anySame)
}
