makeonecross <-
function(x1,x2,Candidates,mutprob=.5){
	n1<-length(unlist(x1))
	n2<-length(unlist(x2))
	n<-min(c(n1,n2))
	x1x2<-union(unlist(x1),unlist(x2))
	cross<-sample(x1x2,n, replace=FALSE)
	randnum<-runif(1)
	if (randnum<mutprob){
		cross[sample(1:n,1)]<-sample(setdiff(Candidates,cross),1)
	}
	return(cross)
}
