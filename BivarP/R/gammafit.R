# require(dfoptim)
gammafit <- function(time, event) {
  gammaLik <- function(x) {
    ie <- which(event==0)
    it <- which(event==1)
    if(length(it)>0) dg <- dgamma(time[it], shape=x[1], scale=x[2], log=TRUE)
    else print("count of noncensored times is 0")
    if(length(ie) > 0) sg <- pgamma(time[ie], shape=x[1], scale=x[2],
      lower.tail=FALSE, log.p=TRUE) else sg = 0
    return(-(sum(dg)+sum(sg)))
  }
  no = length(time)
  s = log(sum(time)/no) - sum(log(time))/no
  sx = (3 - s + sqrt((s-3)*(s-3) + 24*s))/(12*s)
  xx = c(sx, max(time)/2)
  return(nmkb(xx, gammaLik, lower=c(0.01, 0.1), upper=c(sx*10, max(time))))
}
