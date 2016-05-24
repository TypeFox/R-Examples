dispersion.all <-
function(x,y,z,delta,trials=c(100,0.25),
	conf=c(0.025,0.975),outfile="Dispersion.pdf")
{
{
	predicted.all<-list(c(1:ncol(y)))
	x.values<-seq(min(x),max(x),delta)
	limits<-as.list(c(1:ncol(y)))
	names(limits)<-colnames(y)
	limits1<-matrix(nrow=(length(x.values)),ncol=5)
	colnames(limits1)<-c("Uppest","Upper","Mean","Lower","Lowest")
	rownames(limits1)<-x.values
	for(i in 1:ncol(y)){
		dispersion(x,y[,i],delta=delta,span=z[i,1],degree=z[i,2],
			trials=trials)->predicted.all[[i]]
		limits1->limits[[i]]
		apply(predicted.all[[i]],1,max)->limits[[i]][,1]
		apply(predicted.all[[i]],1,min)->limits[[i]][,5]
		apply(predicted.all[[i]],1,mean)->limits[[i]][,3]
		apply(predicted.all[[i]],1,quantile,
			conf[1])->limits[[i]][,4]
		apply(predicted.all[[i]],1,quantile,
			conf[2])->limits[[i]][,2]
		}
	loess.res<-as.list(c(1:ncol(y)))
	for(i in 1:ncol(y)){
		loess(y[,i]~x,span=z[i,1],degree=z[i,2])->loess.res[[i]]
		}
	coef.det<-matrix(ncol=3,nrow=ncol(y))
	colnames(coef.det)<-c("SS regression","SS total","R2")
	rownames(coef.det)<-colnames(y)
	for(i in 1:ncol(y)){
		sum((predict(loess.res[[i]],x)-mean(y[,
			i]))^2)->coef.det[i,1]
		sum((y[,i]-mean(y[,i]))^2)->coef.det[i,2]
		}
	coef.det[,3]<-coef.det[,1]/coef.det[,2]
	coef.det<-round(coef.det,digits=2)
	pdf(outfile)
	layout(cbind(c(1,4,7),c(2,5,8),c(3,6,9)))
	for(i in 1:length(predicted.all)){
		plot(x.values,predicted.all[[i]][,1],type="l",xlab="",
			ylab="",main=colnames(y)[i],ylim=c(0,
			max(limits[[i]])))
			polygon(as.vector(c(x.values,rev(x.values))),
				as.vector(c(limits[[i]][,5],rev(limits[[i]][,
				1]))),border=NA,density=-1,col="gray")
			lines(x.values,limits[[i]][,3],type="l")
			lines(x.values,limits[[i]][,2],type="l",lty=2,lwd=0.5)
			lines(x.values,limits[[i]][,4],type="l",lty=2,lwd=0.5)
			points(x,y[,i],pch=18)
			points(x[which(y[,i]>max(limits[[i]]))],
				rep(max(limits[[i]]),length(which(y[,
				i]>max(limits[[i]])))),pch=3,col="red")
			}
	for(i in 1:length(loess.res)){
		qqnorm(loess.res[[i]]$residuals,main=colnames(y)[i])
		qqline(loess.res[[i]]$residuals,lty=2)
		abline(0,1)
		}
	dev.off()
}
results<-list(limits,coef.det)
names(results)<-c("limits","coef.det")
return(results)
}

