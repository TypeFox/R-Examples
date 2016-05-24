compute.RQ<-function(RAN)
{
r.star<-apply(RAN,1,mean)
q.star<-rep(NA,nrow(RAN))
for (i in 1:nrow(RAN))
q.star[i]<-sum((RAN[i,]-r.star[i])^2)
return(cbind(r.star,q.star))
}
