is.data <- function(data){
  if(!is.data.frame(data)){stop("The data set you entered is not a data frame!")}
  for(i in 2:ncol(data)){
    if(!is.numeric(data[[i]])){stop("One or more columns in your data set are not numeric!")}
  }
}