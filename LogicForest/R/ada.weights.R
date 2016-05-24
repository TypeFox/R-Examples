ada.weights <-
function(fit, cwts, Xs, resp)
{
 n<-nrow(Xs)
 prd<-predict(fit, newbin=Xs)
 miss<-ifelse(prd!=resp, 1, 0)
 actual.miss<-sum(miss)/n
 err.num<-cwts*miss
 err.bar<-sum(err.num)/sum(cwts)
 if (err.bar==0) {alpha<-length(cwts)}
 if (err.bar==1) {alpha<--length(cwts)}
 if (err.bar!=0 & err.bar!=1) {alpha<-log((1-err.bar)/err.bar)}
 nwts<-cwts*exp(alpha*miss)
 ans<-list(nwts=nwts, alpha=alpha, actual.miss=actual.miss, wted.miss=err.bar)
}
