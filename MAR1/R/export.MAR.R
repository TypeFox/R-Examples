export.MAR<-function(model.out,export=TRUE){
if(export==T) export<-paste("MAR.results",Sys.time())
if(file.exists(export)) export<-paste(export,Sys.time())

dir.create("MAR.results")

temp<-model.out$variables.selected
write.table(temp,"MAR.results/variables.selected.csv",col.names=F,sep=",")

temp<-model.out$restrictions.set
temp<-data.frame(rownames(temp),stack(as.data.frame(temp)))
temp<-temp[,c(1,3,2)]
names(temp)<-c("response","predictor","restriction")
write.table(temp,"MAR.results/restrictions.set.csv",row.names=F,sep=",")

temp<-as.matrix(model.out$search.type)
colnames(temp)<-"search.type"
write.table(temp,"MAR.results/search.type.csv",row.names=F,sep=",")

temp<-as.matrix(model.out$search.time)
colnames(temp)<-"search.time"
write.table(temp,"MAR.results/search.time.csv",row.names=F,sep=",")


#######
write.results<-function(r,r.name){
r.A<-r$A
write.table(r.A,paste("MAR.results/",r.name,".A.csv",sep=""),col.names=F,sep=",")

r.B<-r$B
r.B<-data.frame(rownames(r.B),stack(as.data.frame(r.B)))
r.B<-r.B[,c(1,3,2)]
names(r.B)<-c("response","predictor","b")
write.table(r.B,paste("MAR.results/",r.name,".B.csv",sep=""),row.names=F,sep=",")

if(length(r$C)>0){
r.C<-r$C
r.C<-data.frame(rownames(r.C),stack(as.data.frame(r.C)))
r.C<-r.C[,c(1,3,2)]
names(r.C)<-c("response","predictor","c")
write.table(r.C,paste("MAR.results/",r.name,".C.csv",sep=""),row.names=F,sep=",")
}

r.lnlk<-r$log.likelihood
names(r.lnlk)<-"log.likelihood"
write.table(t(r.lnlk),paste("MAR.results/",r.name,".log.likelihood.csv",sep=""),row.names=F,sep=",")

r.aic<-r$AIC
names(r.aic)<-"AIC"
write.table(t(r.aic),paste("MAR.results/",r.name,".AIC.csv",sep=""),row.names=F,sep=",")

r.bic<-r$BIC
names(r.bic)<-"BIC"
write.table(t(r.bic),paste("MAR.results/",r.name,".BIC.csv",sep=""),row.names=F,sep=",")

r.r2<-r$R2.values
r.r2<-data.frame("variate"=rownames(r.r2),r.r2)
write.table(r.r2,paste("MAR.results/",r.name,".R2.values.csv",sep=""),row.names=F,sep=",")

r.sdist<-r$stationary.distribution
r.sdist.name<-paste(r.name,"stationary.distribution",sep=".")

temp<-r.sdist$mean
temp<-data.frame(variate=rownames(temp),temp)
write.table(temp,paste("MAR.results/",r.sdist.name,".mean.csv",sep=""),row.names=F,sep=",")

temp<-r.sdist$covariance
temp<-data.frame(rownames(temp),stack(as.data.frame(temp)))
temp<-temp[,c(1,3,2)]
names(temp)<-c("response","predictor","V")
write.table(temp,paste("MAR.results/",r.sdist.name,".covariance.csv",sep=""),row.names=F,sep=",")

r.perr<-r$process.errors
r.perr.name<-paste(r.name,"process.errors",sep=".")

temp<-r.perr$residuals
temp<-data.frame(time.step=1:nrow(temp),temp)
write.table(temp,paste("MAR.results/",r.perr.name,".residuals.csv",sep=""),row.names=F,sep=",")

temp<-r.perr$covariance
temp<-data.frame(rownames(temp),stack(as.data.frame(temp)))
temp<-temp[,c(1,3,2)]
names(temp)<-c("response","predictor","sigma")
write.table(temp,paste("MAR.results/",r.perr.name,".covariance.csv",sep=""),row.names=F,sep=",")

temp<-r.perr$corrmatrix
temp<-data.frame(rownames(temp),stack(as.data.frame(temp)))
temp<-temp[,c(1,3,2)]
names(temp)<-c("response","predictor","correlation")
write.table(temp,paste("MAR.results/",r.perr.name,".corrmatrix.csv",sep=""),row.names=F,sep=",")

r.stab<-r$stability
r.stab.name<-paste(r.name,"stability",sep=".")

temp<-matrix(r.stab$resilience$eigB)
colnames(temp)<-"eigB"
write.table(temp,paste("MAR.results/",r.stab.name,".resilience.eigB.csv",sep=""),row.names=F,sep=",")

temp<-r.stab$resilience$detB
names(temp)<-"detB"
write.table(t(temp),paste("MAR.results/",r.stab.name,".resilience.detB.csv",sep=""),row.names=F,sep=",")

temp<-r.stab$resilience$maxeigB
names(temp)<-"maxeigB"
write.table(t(temp),paste("MAR.results/",r.stab.name,".resilience.maxeigB.csv",sep=""),row.names=F,sep=",")

temp<-r.stab$resilience$maxeigkrB
names(temp)<-"maxeigkrB"
write.table(t(temp),paste("MAR.results/",r.stab.name,".resilience.maxeigkrB.csv",sep=""),row.names=F,sep=",")

temp<-r.stab$reactivity$sigma.over.Vinf
names(temp)<-"sigma.over.Vinf"
write.table(t(temp),paste("MAR.results/",r.stab.name,".reactivity.sigma.over.Vinf.csv",sep=""),row.names=F,sep=",")

temp<-r.stab$reactivity$maxeigBxB
names(temp)<-"maxeigBxB"
write.table(t(temp),paste("MAR.results/",r.stab.name,".reactivity.maxeigBxB.csv",sep=""),row.names=F,sep=",")
}

write.results(r=model.out$bestfit,r.name="bestfit")

if(!is.null(model.out$bootstrap)) {
	write.results(r=model.out$bootstrap,r.name="bootstrap")
	
r<-model.out$bootstrap$limits

temp<-data.frame(variate=rownames(r$bootA),lower=r$lowerA,upper=r$upperA)
write.table(temp,"MAR.results/bootstrap.limits.A.csv",row.names=F,sep=",")

temp<-data.frame(rownames(r$bootB),stack(as.data.frame(r$bootB)))
names(temp)<-c("response","b","predictor")
temp$lower<-stack(as.data.frame(r$lowerB))[,1]
temp$upper<-stack(as.data.frame(r$upperB))[,1]
temp<-temp[,c(1,3,4,5)]
write.table(temp,"MAR.results/bootstrap.limits.B.csv",row.names=F,sep=",")

if(length(r$bootC)>0){
temp<-data.frame(rownames(r$bootC),stack(as.data.frame(r$bootC)))
names(temp)<-c("response","c","predictor")
temp$lower<-stack(as.data.frame(r$lowerC))[,1]
temp$upper<-stack(as.data.frame(r$upperC))[,1]
temp<-temp[,c(1,3,4,5)]
write.table(temp,"MAR.results/bootstrap.limits.C.csv",row.names=F,sep=",")
}

r<-model.out$bootstrap$limits$stationary.distribution

temp<-data.frame(variate=rownames(r$boot.mean),lower=r$lower.mean,upper=r$upper.mean)
write.table(temp,"MAR.results/bootstrap.limits.stationary.distribution.mean.csv",row.names=F,sep=",")

temp<-data.frame(rownames(r$boot.covariance),stack(as.data.frame(r$boot.covariance)))
names(temp)<-c("response","mean","predictor")
temp$lower<-stack(as.data.frame(r$lower.covariance))[,1]
temp$upper<-stack(as.data.frame(r$upper.covariance))[,1]
temp<-temp[,c(1,3,4,5)]
write.table(temp,"MAR.results/bootstrap.limits.stationary.distribution.covariance.csv",row.names=F,sep=",")

r<-model.out$bootstrap$limits$process.errors
temp$lower<-as.vector(r$lower.sigma)
temp$upper<-as.vector(r$upper.sigma)
write.table(temp,"MAR.results/bootstrap.limits.process.errors.csv",row.names=F,sep=",")

	}

if(length(model.out$top.bestfit)>0){
	r<-model.out$top.bestfit
	temp<-data.frame(response=rep(rownames(r),ncol(r)),predictor=rep(colnames(r),each=nrow(r)))
	for(i in 1:dim(r)[3]) temp<-cbind(temp,stack(as.data.frame(r[,,i]))[,1])
	names(temp)[-c(1:2)]<-paste("AIC_",dimnames(r)[[3]],sep="")
	write.table(temp,"MAR.results/top.bestfit.csv",row.names=F,sep=",")
	}

#######

invisible(file.rename("MAR.results",gsub(":",".",export)))
}

