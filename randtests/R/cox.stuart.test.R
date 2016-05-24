##
##  Cox Stuart Test
##
cox.stuart.test <- function(x, alternative="two.sided"){
  # Performs the Cox Stuart trend test.
  #
  # Args:
  #   x: data values.
  #
  # Returns:
  #   statistic: the value of the statistic test.
  #   n: the sample size, after the remotion of observations with same value.   
  #   p.value: p-value.
  #
  dname <- deparse(substitute(x))
  if (alternative == "t"){alternative <- "two.sided"} 
  if (alternative == "l"){alternative <- "left.sided"}
  if (alternative == "r"){alternative <- "right.sided"}    
  if (alternative != "two.sided" & alternative != "left.sided" & alternative != "right.sided")
  {stop("must give a valid alternative")}    
  stopifnot(is.numeric(x))  
  n0 <- length(x)
  if (n0 < 2){stop("sample size must be greater than 1")} 
  n0<-round(length(x)) %% 2
  if (n0 == 1) {
    remove <- (length(x)+1)/2
    x <- x[ -remove ] 
  } 
  half <- length(x)/2
  x1 <- x[1:half]
  x2 <- x[(half+1):(length(x))]
  n<-sum((x2-x1)!=0)
  t<-sum(x1<x2)
  if (alternative=="left.sided") {
    p.value<-pbinom(t,n,0.5)
    alternative<-"decreasing trend"
  }
  if (alternative=="right.sided"){p.value<- 1-pbinom(t-1,n,0.5);alternative<-"increasing trend"}
  if (alternative=="two.sided"){p.value<-2*(min(1-pbinom(t-1,n,0.5),pbinom(t,n,0.5)));alternative<-"non randomness"}
  #pvalue<-min(pvalue,1)
  #names(n)="n"
  rval<-list(statistic=c(statistic=t), alternative=alternative, p.value=p.value,
            method="Cox Stuart test", parameter=c(n=n), data.name=dname)
  class(rval)<-"htest"
  return(rval)
}
  
  
  
  
 