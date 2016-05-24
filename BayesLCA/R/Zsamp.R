Zsamp <-
function(prob, counts)
{
	G<-ncol(prob);z<-array(0,dim(prob), dimnames=list(NULL,1:G)) 
	for(n in 1:nrow(prob)){
		tab<-table(sample(1:G, counts[n], replace=TRUE, prob[n,]));z[n,names(tab)]<-tab}
		z
	}
