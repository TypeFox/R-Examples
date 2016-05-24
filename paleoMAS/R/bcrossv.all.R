bcrossv.all <-
function(x,y,interval=c(0.15,1,0.05),
	trials=c(10,0.25),target=c("rse","rmse"))
{
{
	crossval<-matrix(nrow=ncol(y),ncol=4)
	colnames(crossval)<-c("alpha","degree","rse","rmse")
	rownames(crossval)<-colnames(y)
	crossval1<-as.list(c(1:ncol(y)))
	if(target=="rse"){
	for(i in 1:ncol(y)){
		bcrossv.l(x,y[,i],interval=interval,
			trials=trials,plot=FALSE)->crossval1[[i]]
		crossval1[[i]][which(crossval1[[i]][,
				3]==max(crossval1[[i]][,3]))[1],1:4]->crossval[i,]
		}
	}
	if(target=="rmse"){
	for(i in 1:ncol(y)){
		bcrossv.l(x,y[,i],interval=interval,
			trials=trials,plot=FALSE)->crossval1[[i]]
		crossval1[[i]][which(crossval1[[i]][,
				4]==max(crossval1[[i]][,4]))[1],1:4]->crossval[i,]
		}
	}
}
return(crossval)
}

