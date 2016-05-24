stat.gen <-
function(num.subj, target.mean, target.sd, data.dec=2, neg.data=TRUE, x.name="Variable X"){
	iteration<-0
	target.mean<-round(target.mean,data.dec+1)
	target.sd<-round(target.sd,data.dec+1)
	dif<-target.mean-target.sd	
		repeat{
			iteration=iteration+1
			cat("working on iteration",iteration,"\r")
			ifelse(neg.data==FALSE & dif<=0.5,x<-scale(rexp(num.subj)),x<-scale(rnorm(num.subj,target.mean,target.sd)))
			ifelse(neg.data==FALSE,x.2<-abs(round(matrix(x[,1]*target.sd+target.mean),data.dec)),x.2<-round(matrix(x[,1]*target.sd+target.mean),data.dec))			
			if(round(apply(x.2,2,mean),data.dec+1)==target.mean&round(apply(x.2,2,sd),data.dec+1)==target.sd|iteration==100000)
			break
			}
	if(iteration==100000)	cat("Maximum Number of Iterations Reached Without A Solution. Try Again!","\b")
		else{
			x.2<-round(matrix(x[,1]*target.sd+target.mean),data.dec)
			colnames(x.2)<-x.name
			write.table(x.2,file=paste(x.name,"MorseGen Results.txt",sep=" "),row.names=FALSE)
			sim.perf<-data.frame(rbind(c(target.mean,round(apply(x.2,2,mean),data.dec+2),target.sd,round(apply(x.2,2,sd),data.dec+2),iteration)))
			colnames(sim.perf)<-c(paste("Target",x.name,"Mean",sep=" "),paste("Sample",x.name,"Mean",sep=" "),paste("Target",x.name,"SD",sep=" "),paste("Sample",x.name,"SD",sep=" "),"Iterations")
			cat("Summary of the Simulation Results (note: raw data has been exported to the file: ")
			cat(x.name)
			cat(" MorseGen Results.txt in your working directory","\b")
			cat(")","\n")
			list(DATA=x.2, PERFORMANCE=sim.perf)
			}
}
