summary.MAR<-function(object,...){

	# msummary<-list()

make.coef.sum<-function(object.which){
	# attach(object.which)
	# on.exit(detach(object.which))
	b<-data.frame(matrix="B",coef=as.vector(object.which$B))
	if(length(object$bestfit$C)==0){
	c<-NULL } else {
	c<-data.frame(matrix="C",coef=as.vector(object.which$C))
	colnames(c)<-colnames(b)}
	
	coef.all<-rbind(b,c)

	coef.sum<-data.frame(matrix=unique(coef.all[,1]))
	matnames<-data.frame(matrix=unique(coef.all[,1]))

	coef.sum$total.coef	<-merge(aggregate(coef~matrix,coef.all,length),matnames,all=T,sort=T)[,2]	
	if(nrow(coef.all[coef.all$coef==0,])>0){
	coef.sum$zeros		<-merge(aggregate(coef~matrix,coef.all[coef.all$coef==0,],length),matnames,all=T,sort=T)[,2]
	}else{coef.sum$zeros<-0}
	if(nrow(coef.all[coef.all$coef!=0,])>0){
	coef.sum$nonzeros		<-merge(aggregate(coef~matrix,coef.all[coef.all$coef!=0,],length),matnames,all=T,sort=T)[,2]
	}else{coef.sum$nonzeros<-0}
	if(nrow(coef.all[coef.all$coef>0,])>0){
	coef.sum$positive		<-merge(aggregate(coef~matrix,coef.all[coef.all$coef>0,],length),matnames,all=T,sort=T)[,2]
	}else{coef.sum$positive<-0}
	if(nrow(coef.all[coef.all$coef<0,])>0){
	coef.sum$negative		<-merge(aggregate(coef~matrix,coef.all[coef.all$coef<0,],length),matnames,all=T,sort=T)[,2]
	}else{coef.sum$negative<-0}

	coef.sum[is.na(coef.sum)]<-0
	coef.sum[,1]<-factor(coef.sum[,1],levels=c("B","C","Total"))
	coef.sum[3,]<-data.frame(matrix="Total",data.frame(t(apply(coef.sum[,-1],2,sum))))
	if(is.na(coef.sum[2,1])) coef.sum[2,1]<-"C"
	
	rownames(coef.sum)<-coef.sum[,1]
	coef.sum<-t(coef.sum[,-1])
	
	coef.sum
	}
	
coef.sum.best<-make.coef.sum(object$bestfit)
ics.best<-c(object$bestfit$AIC,object$bestfit$BIC)
names(ics.best)<-c("AIC","BIC")
r2.best<-summary(object$bestfit$R2.values)
	r2.best<-apply(substr(r2.best,nchar("1st Qu.: "),max(nchar(r2.best))),c(1,2),as.numeric)
	colnames(r2.best)<-gsub(" ","",colnames(r2.best))
	rownames(r2.best)<-substr(summary(object$bestfit$R2.values)[,1],1,nchar("1st Qu"))
	rownames(r2.best)<-gsub("[.]","",rownames(r2.best))
stab.temp<-object$bestfit$stability
stab.best<-data.frame( Attribute=c(rep("resilience",3),rep("reactivity",2)),
	Metric=c(names(stab.temp[[1]])[-1],names(stab.temp[[2]])),
	Value=c(unlist(stab.temp[[1]][-1]),unlist(stab.temp[[2]])) )
	rownames(stab.best)<-NULL

if(!is.null(object$bootstrap)) {
	coef.sum.boot<-make.coef.sum(object$bootstrap)
	ics.boot<-c(object$bootstrap$AIC,object$bootstrap$BIC);names(ics.boot)<-c("AIC","BIC")
	r2.boot<-summary(object$bootstrap$R2.values)
		r2.boot<-apply(substr(r2.boot,nchar("1st Qu.: "),max(nchar(r2.boot))),c(1,2),as.numeric)
		colnames(r2.boot)<-colnames(r2.best)
		rownames(r2.boot)<-rownames(r2.best)
	stab.temp<-object$bootstrap$stability
	stab.boot<-data.frame( Attribute=c(rep("resilience",3),rep("reactivity",2)),
		Metric=c(names(stab.temp[[1]])[-1],names(stab.temp[[2]])),
		Value=c(unlist(stab.temp[[1]][-1]),unlist(stab.temp[[2]])) )
		rownames(stab.boot)<-NULL

	} else {
		coef.sum.boot<-NULL
		ics.boot<-NULL
		r2.boot<-NULL
		stab.boot<-NULL
		}
	

	msummary<-list(coef.summary.best=coef.sum.best,coef.summary.boot=coef.sum.boot,
		ICs.best=ics.best,ICs.boot=ics.boot,
		R2.best=r2.best,R2.boot=r2.boot,
		stability.best=stab.best,stability.boot=stab.boot)
	
	class(msummary)<-"MARsummary"
	msummary

}
