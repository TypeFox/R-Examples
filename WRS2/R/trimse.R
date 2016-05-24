trimse <- function(x, tr = 0.2, na.rm = FALSE){
#
#  Estimate the standard error of the gamma trimmed mean
#  The default amount of trimming is tr=.2.
#
if(na.rm)x<-x[!is.na(x)]
trimse<-sqrt(winvar(x,tr))/((1-2*tr)*sqrt(length(x)))
trimse
}

