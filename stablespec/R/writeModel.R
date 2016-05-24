writeModel <- function(theModel, numVar, longitudinal) {

  if (longitudinal) {
    numVar <- numVar * 2
  }

  indeks <- 0
  forVariance <- modelString <- NULL

  for (j in 1:numVar) {
      counter <- j
      jExists <- FALSE

    for (k in 1:numVar) {
        if (j != k && counter == j) {
          #rowwise; row by row
          if (theModel[k, j] == 1){
          jExists <- TRUE
          indeks <- indeks + 1
          modelString <- paste(modelString, 'var', j, ' = par',
                               indeks , '*var', k, ' ', sep="")
          counter <- counter + 1
        }
      } else if (j != k && counter != j) {
        if (theModel[k,j] == 1){
            indeks <- indeks + 1
            modelString <- paste(modelString, '+ par', indeks,
                                 '*var', k, ' ', sep="")
            }
      }
    }

    if (jExists) {
      modelString <- paste(modelString, '\n')
    } else {
      forVariance <- c(forVariance, j)
    }
  }

  #for variances
  for (i in 1:length(forVariance)) {
      modelString <- paste(modelString, 'V(var', forVariance[i],
                           ') = variance', i, '\n', sep="")
      }
  return(modelString)
}
