mWindow <-
function (data, dimens, alpha = 1, windowSize = 2) {

  if (!is.data.frame(data)) {
    stop ("class(data) != 'data.frame'")
  }
  else if (is.null(colnames(data))) {
    stop ("Please add column names to data")
  }
  
  toKeep <- complete.cases(data) 
  data <- data[toKeep,,drop = F]

  check_args_mWindow (data, dimens, alpha, windowSize)

  return(mWindow_main (data, dimens, alpha, windowSize))

}
