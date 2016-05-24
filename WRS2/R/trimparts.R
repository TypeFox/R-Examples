trimparts <-
function(x,tr=.2){
#
#  Compute the trimmed mean, effective sample size, and squared standard error.
#  The default amount of trimming is tr=.2.
#
#  This function is used by other functions described in chapter 6.
#
tm<-mean(x,tr)
h1<-length(x)-2*floor(tr*length(x))
sqse<-(length(x)-1)*winvar(x,tr)/(h1*(h1-1))
trimparts<-c(tm,sqse)
trimparts
}
