cohenW<-
function (X) 
{
n<-sum(X)
W<-sqrt(chisq.test(X)$statistic/n)
as.numeric(W)
}