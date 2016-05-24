bcrossv.l <-
function(x,y,interval=c(0.15,1,0.05),
	trials=c(10,0.25),plot=TRUE)
{
library(lattice)
{
	alphas1<-seq(interval[1],interval[2],interval[3])
	alphas<-rep(alphas1,2)
	degrees<-rep(c(1,2),each=length(alphas1))
	param<-cbind(alphas,degrees)
	crossval<-matrix(nrow=length(alphas),ncol=4)
	colnames(crossval)<-c("alpha","degree","rse","rmse")
	crossval[,1]<-alphas
	crossval[,2]<-degrees
	for(i in 1:nrow(crossval)){
		bcrossv.l1(x,y,trials=trials,span=param[i,1],
			degree=param[i,2],
			plot=FALSE)[2,1]->crossval[i,3]
		bcrossv.l1(x,y,trials=trials,span=param[i,1],
			degree=param[i,2],
			plot=FALSE)[3,1]->crossval[i,4]
		}
	crossval<-round(crossval,2)	
	if(plot==TRUE){
		x<-crossval
		dataplot<-as.data.frame(matrix(nrow=nrow(x)*2,ncol=4))
		dataplot[,1:2]<-rbind(x[,c(1,3)],x[,c(1,4)])
		a<-rep(x[,2],2)
		dataplot[,3]<-ifelse(a==1,"degree 1","degree 2")
		dataplot[,4]<-rep(c("rse","rmse"),each=nrow(x))
		colnames(dataplot)<-c("alpha","value","degree",
			"statistic")
		a<-xyplot(value~alpha|degree*statistic,data=dataplot,
			type="l")
		print(a)
		}
}
return(crossval)
}

