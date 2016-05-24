library(rCUR)

A<-matrix(runif(900),nrow=30);

nr=integer(); nc=integer()
for (i in 1:10000) {
  CUR_A<-CUR(A,k=5,c=2,r=2);
  nr[i]=length(CUR_A@R.index)
  nc[i]=length(CUR_A@C.index)
  if(abs(sum(CUR_A@C.leverage.score)-1)>.Machine$double.eps*10) stop("leverage scores does not sum to 1")
  if(abs(sum(CUR_A@R.leverage.score)-1)>.Machine$double.eps*10) stop("leverage scores does not sum to 1")
}
if(mean(nc)>3 || mean(nr)>3) stop("expected number of selected rows/columns too large")


