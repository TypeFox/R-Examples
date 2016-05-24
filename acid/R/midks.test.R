midks.test <-
function (x, y, ...,w=NULL,pmt=NULL){
  DNAME <- deparse(substitute(x))
  x <- x[!is.na(x)]
  n <- length(x)
  if(is.null(w)) w<-rep(1,n) # w if not defined set to 1 for every observation
  if(is.null(pmt)) pmt<-0 # set to 0 if not defined
  if (n < 1L) stop("not enough 'x' data")
  if (!is.function(y)) stop("'y' must be numeric or a function or a string naming a valid function")
  METHOD <- "One-sample Kolmogorov-Smirnov test, two.sided"
  w<-w/sum(w) # coerce weights to sum to one
  emp.cdf<- c(0,cumsum(w[order(x)]))
  diffpm<-abs(sum(w[x==pmt])-cdf.mix.dag(pmt,...)) # difference for point mass
  x.cont<-x[x>pmt] # x is great than point mass threshold - otherwise problems with upper and lower bound
  w.cont<-w[x>pmt]
  emp.cdf<- c(sum(w[x==pmt]),cumsum(w[order(x.cont)])+sum(w[x==pmt])) # empirical cdf above threshold0
  n.cont<-length(x.cont)
  if(n.cont>=1){
    diff1<-abs(y(sort(x.cont), ...) - emp.cdf[-(n.cont+1)]) # difference between para cdf and empirical cdf, lower, point mass excluded
    diff2<-abs(emp.cdf[-1]-y(sort(x.cont), ...)) # difference between para cdf and empirical cdf, higher, point mass excluded  
  } else diff1<-diff2<-0
  STATISTIC<-max(c(diffpm,diff1,diff2))  
  list(statistic = STATISTIC,method = METHOD,diffpm=diffpm,diff1=diff1,diff2=diff2)
}
