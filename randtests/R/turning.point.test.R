## x<-c(1,3,3,4,6,8,9,3,5,6,7,8,9,10,11,12, 13, 14,15)
##  turning point code
##
turning.point.test <- function(x, alternative="two.sided"){
  stopifnot(is.numeric(x))
  dname <- deparse(substitute(x))
   # Delete consecutive repeated values
  if (min(abs(diff(x)))==0){
    d<-which(diff(x) %in% 0)
    x<-x[-d]
  }
    # Main code  
  n <- length(x)
  mu <- 2*(n-2)/3
  var <- (16*n-29)/90
  
  tp <- sign(apply(embed(x,3), 1, diff))
  tp <- tp[1,]*tp[2,]
  test.sum <- -sum(tp[tp<0])
  test <- (test.sum-mu)/sqrt(var)
  # p-value
  pv0 <- pnorm(test)
  if (alternative=="two.sided"){pv <- 2*min(pv0,1-pv0); alternative<-"non randomness"}
  if (alternative=="left.sided"){pv <- pv0; alternative<-"positive serial correlation"}
  if (alternative=="right.sided"){pv <- 1-pv0; alternative<-"negative serial correlation"}
  # output
  names(n)="n"
  rval <- list(statistic = c(statistic=test), tp=test.sum, p.value = pv, method = "Turning Point Test", 
               data.name = dname, parameter=n, n=n, alternative=alternative) 
  class(rval) <- "htest"
  return(rval)
}