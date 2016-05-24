is.thres.error <- function(x){ # test for threshold.TotalError parameter
  if(!is.numeric(x)){stop("The threshold.TotalError you entered is not numeric!")}
  if(!is.na(x[2])){stop("The threshold.TotalError you entered is a vector!")}
  if(!(x>0)){stop("The threshold.TotalError you entered is zero or less!")}
}