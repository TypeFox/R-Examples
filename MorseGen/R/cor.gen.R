cor.gen <-
function(num.subj, x.mean, x.sd, y.mean, y.sd, correlation, data.dec=2, x.name="Variable X", y.name="Variable Y"){
	x<-rnorm(num.subj)
	y<-rnorm(num.subj)
	xy.res<-resid(lm(y~x))
	x.z<-scale(x)
	res.z<-scale(xy.res)
	z<-correlation*x.z+sqrt(1-correlation^2)*res.z
	Variable.X<-x.mean+x.sd*x.z
	Variable.Y<-y.mean+y.sd*z
	data<-round(data.frame(Variable.X,Variable.Y),data.dec)
	colnames(data)<-c(x.name,y.name)
	sim.perf<-data.frame(rbind(c(apply(Variable.X,2,mean),apply(Variable.X,2,sd),apply(Variable.Y,2,mean),apply(Variable.Y,2,sd),cor(Variable.X,Variable.Y))))
	colnames(sim.perf)<-c(paste("Sample",x.name,"Mean",sep=" "),paste("Sample",x.name,"SD",sep=" "),paste("Sample",y.name,"Mean",sep=" "),paste("Sample",y.name,"SD",sep=" "),"  Sample Correlation")
	targets<-data.frame(rbind(c(x.mean,x.sd,y.mean,y.sd,correlation)))
	colnames(targets)<-c(paste("Target",x.name,"Mean",sep=" "),paste("Target",x.name,"SD",sep=" "),paste("Target",y.name,"Mean",sep=" "),paste("Target",y.name,"SD",sep=" "),"  Target Correlation")
	cat("Summary of the Simulation Results (note: raw data has been exported to the file: ")
	cat(x.name)
	cat(" ")
	cat(y.name)
	cat(" Correlation")
	cat(" MorseGen Results.txt in your working directory","\b")
	cat(")","\n")
	write.table(data,file=paste(x.name,y.name,"Correlation MorseGen Results.txt",sep=" "),row.names=FALSE,sep="\t")
	list(DATA=data, TARGETS=targets, RESULTS=sim.perf)
}
