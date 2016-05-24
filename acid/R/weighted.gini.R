weighted.gini <-
function(x,w=NULL){
  if(is.null(w)) w<-rep(1,length(x))
  x.sort<- sort(x)
  x.order<- order(x)
  x<-x.sort
  n<- length(x)
  w<-w[x.order]/sum(w)
  w.c <- cumsum(w)
  xw.c <- cumsum(w*x)
  xw.c <- xw.c/xw.c[n] #coercing such cumulative distr with max=1
  Gini<-t(xw.c[-1])%*%w.c[-n] - t(xw.c[-n])%*%w.c[-1]
  Ginistar <- Gini * n/(n - 1) # bias corrected
  Ginistarw<- Gini * 1/(1-sum(w^2)) # using bias correction from cov.wt
  #  print("Warning: weighting is not properly accounted for in the sample adjustment of bcGini!!")
  list(Gini=Gini,"bcGini"=Ginistar,"bcwGini"=Ginistarw)
}
