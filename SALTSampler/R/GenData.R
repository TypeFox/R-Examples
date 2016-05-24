GenData <- function(center, n, size) {
  #Generates a synthetic data set representing multiple draws from a
  #multinomial distribution with user-specified parameters. A matrix
  #n rows corresponding to each draw is outputted where the entry in 
  #the ith column and the jth row gives the number of the items that 
  #were in the ith bin on the jth trial
  
  #Check inputs
  if (!exists("center")) {
    stop("center is not defined")
  }
  if (!is.numeric(center)){
    stop("center is not numeric")
  }
  if (!is.vector(center)){
    stop("center is not a vector")
  }
  
  if (!exists("n")) {
    stop("n is not defined")
  }  
  if (!is.vector(n)){
    stop("n is not a vector")
  }
  if (n%%1 != 0){
    stop("n is not an integer")
  }
  if (length(n) != 1){
    stop("n is not of length 1")
  }
  
  if (!exists("size")) {
    stop("n is not defined")
  }  
  if (!is.vector(size)){
    stop("n is not a vector")
  }
  if (size%%1 != 0){
    stop("n is not an integer")
  }
  if (length(size) != 1){
    stop("n is not of length 1")
  }
  
  #Rescale center to proportion
  center <- center/sum(center)
  
  #Sample and reformat multinomial data for use in sampler
  data <- t(rmultinom(n, size, center))
  
  #Return data
  return(data)
}
