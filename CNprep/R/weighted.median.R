weighted.median <-
function(v,weights){
	weights<-weights[order(v)]
	v<-sort(v)
	sw<-sum(weights)
	return(v[which.min(abs(cumsum(weights)-0.5*sw))])
}
