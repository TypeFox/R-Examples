sunterpi <-
function(x,n)
{
xo<-sort(x,decreasing=TRUE)
pik<-n*xo/sum(xo)
if(pik[1]>1) stop("There are some units with inclusion probability >1")
N<-length(xo)
t<-rev(cumsum(rev(xo)))
#xbar<-t/(N:1)
kk0<-n*xo/t
k0<-which(kk0>=1,arr.ind = FALSE)[1]
kstar<-min(k0,N-n+1)
pik[kstar:N]<-n*t[kstar]/(N-kstar+1)/sum(x)
ord<-N+1-rank(x,ties.method="random")
pik<-pik[ord]
pik
}

