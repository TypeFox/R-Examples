error.l <-
function(x,y,z,trials=c(100,0.25))
{
{
	crossval<-matrix(nrow=ncol(y),ncol=6)
	rownames(crossval)<-colnames(y)
	colnames(crossval)<-c("rse","rmse","range","mean","rmse/range",
		"rmse/mean")
	for(i in 1:ncol(y)){
		bcrossv.l1(x,y[,i],trials=trials,span=z[i,1],degree=z[i,2],
			plot=FALSE)[2:3,]->crossval[i,1:2]
		}
	crossval[,3]<-apply(y,2,range)[2,]-apply(y,2,range)[1,]
	crossval[,4]<-apply(y,2,mean)
	crossval[,5]<-crossval[,2]/crossval[,3]
	crossval[,6]<-crossval[,2]/crossval[,4]
	crossval<-round(crossval,2)
}
return(crossval)
}

