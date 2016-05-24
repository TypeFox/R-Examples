getwm <-
function(orders,skeleton){
	d<-ncol(orders)
	s<-nrow(orders)
	alpha<-matrix(0,nrow=s,ncol=d)
	for(j in 1:s){
	alpha[j,]<-skeleton[order(orders[j,])]
}
alpha
}
