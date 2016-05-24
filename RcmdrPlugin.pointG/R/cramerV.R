cramerV<-
function (X) 
{
nr<-nrow(X)
nl<-ncol(X)
n<-sum(X)
V<-sqrt(chisq.test(X)$statistic/(n*min(nr-1,nl-1)))
as.numeric(V)
}