is.sampleLeng <- function(x) { # test for sampleLength
  if(!is.numeric(x)){stop("The sampleLength you entered is not numeric!")}
  if(!is.na(x[2])){stop("The sampleLength you entered is a vector!")}
  if(!(x>0)){stop("The sampleLength you entered is zero or less!")}
  if(!(x<1)){stop("The sampleLength you entered is one or more!")}
}