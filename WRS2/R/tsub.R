tsub<-function(isub,x,y,tr){
#
#  Compute test statistic for trimmed means
#  when comparing dependent groups.
#  By default, 20% trimmed means are used.
#  isub is a vector of length n,
#  a bootstrap sample from the sequence of integers
#  1, 2, 3, ..., n
#
#  This function is used by ydbt
#
tsub<-yuend(x[isub],y[isub],tr=tr)$test
tsub
}
