yuend <- function(x, y, tr = 0.2){
#
#  Compare the trimmed means of two dependent random variables
#  using the data in x and y.
#  The default amount of trimming is 20%
#
#  Any pair with a missing value is eliminated
#  The function rm2miss allows missing values.
#
#  A confidence interval for the trimmed mean of x minus the
#  the trimmed mean of y is computed and returned in yuend$ci.
#  The significance level is returned in yuend$p.value
#
#  This function uses winvar from chapter 2.
#
  cl <- match.call()
  alpha=.05  
  if(length(x)!=length(y))stop("The number of observations must be equal")
  m<-cbind(x,y)
  m<-elimna(m)
  x<-m[,1]
  y<-m[,2]
  h1<-length(x)-2*floor(tr*length(x))
  q1<-(length(x)-1)*winvar(x,tr)
  q2<-(length(y)-1)*winvar(y,tr)
  q3<-(length(x)-1)*wincor(x,y,tr)$cov
  df<-h1-1
  se<-sqrt((q1+q2-2*q3)/(h1*(h1-1)))
  crit<-qt(1-alpha/2,df)
  dif<-mean(x,tr)-mean(y,tr)
  low<-dif-crit*se
  up<-dif+crit*se
  test<-dif/se
  yuend<-2*(1-pt(abs(test),df))
  
  result <- list(test = test, conf.int = c(low,up), se = se, p.value = yuend, df = df, diff = dif, call = cl)
  class(result) <- "yuen"
  result
}
