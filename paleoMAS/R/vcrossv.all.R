vcrossv.all <-
function(x,f,to,nsimulat,funct,ntrials,plot=TRUE)
{
{
	accuracy1<-matrix(ncol=to,nrow=ntrials)
	x1<-as.list(c(1:ntrials))
	f1<-as.list(c(1:ntrials))
	a<-vcrossv.da(x,f,fold=1,nsimulat=nsimulat,
		funct=funct)$accuracy
	accuracy1[,1]<-rep(a,ntrials)
	for(j in 2:ncol(accuracy1)){
		for(i in 1:nrow(accuracy1)){
			x1[[i]]<-x
			f1[[i]]<-f
			vcrossv.da(x1[[i]],f1[[i]],fold=j,nsimulat=nsimulat,
				funct=funct)$accuracy->accuracy1[i,j]
			}
		}
	accuracy<-matrix(nrow=ncol(accuracy1),ncol=4)
	colnames(accuracy)<-c("fold","mean accuracy",
		"lower (0.025)","upper (0.975)")
	accuracy[,1]<-c(1:to)
	accuracy[,2]<-apply(accuracy1,2,mean)
	accuracy[,3]<-apply(accuracy1,2,quantile,probs=0.025)
	accuracy[,4]<-apply(accuracy1,2,quantile,probs=0.975)
	if(plot==TRUE){
		plot.data<-matrix(nrow=to,ncol=4)
		plot.data[,1]<-c(1:to)
		plot.data[,2]<-apply(accuracy1,2,min)
		plot.data[,3]<-apply(accuracy1,2,mean)
		plot.data[,4]<-apply(accuracy1,2,max)
		plot(accuracy[,c(1,2)],ylim=c(0,max(accuracy1)),type="n")
		arrows(plot.data[-1,1],plot.data[-1,3],plot.data[-1,
			1],plot.data[-1,2],length=0.05,angle=90)
		arrows(plot.data[-1,1],plot.data[-1,3],plot.data[-1,
			1],plot.data[-1,4],length=0.05,angle=90)
		lines(accuracy[,c(1,2)])
		lines(accuracy[,c(1,3)],lty=2,col="gray")
		lines(accuracy[,c(1,4)],lty=2,col="gray")
		}
}
return(accuracy)
}

