
myWdataclean2<-function (z, d, wt = rep(1, length(z))) 
{
#sort z,d,wt (ascending) wrt z, using -d to order ties 
    niceorder <- order(z, -d)
    sortedz <- z[niceorder]
    sortedd <- d[niceorder]
    sortedw <- wt[niceorder]
#store length of sortedd in n
    n <- length(sortedd)
#y1 checks for jumps in sortedz using offsets of sortedz
    y1 <- sortedz[-1] != sortedz[-n]
#y2 checks for "jumps" in sortedd 
    y2 <- sortedd[-1] != sortedd[-n]
#y checks for jumps (in sortedz or sortedd) using y1 and y2
    y <- y1 | y2
#ind stores jump indices (final index will be  n)
    ind <- c(which(y | is.na(y)), n)
#csum is cumulative sum of the weights
    csumw <- cumsum(sortedw)
#value contains the (unique) obs in sortedz
#dd contains the status of obs in sortedz
#weight has the weights of the obs in sortedz
    list(value = sortedz[ind], dd = sortedd[ind], weight = diff(c(0,csumw[ind])))
}