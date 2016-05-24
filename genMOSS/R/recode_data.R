recode_data <-
function (data, dimens, alpha = 1) {

  if (!is.data.frame(data)) {
    stop ("class(data) != 'data.frame'")
  }
  else if (is.null(colnames(data))) {
    stop ("Please add column names to data")
  }
  
  toKeep <- complete.cases(data) 
  data <- data[toKeep,,drop = F]

  check_args_recode_data (data, dimens, alpha) 

  return(recode_data_main(data, dimens, alpha))

}
