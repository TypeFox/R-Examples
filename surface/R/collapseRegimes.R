collapseRegimes <-
function(otree, odata, oldshifts, oldaic, oldfit, aic_threshold=0, only_best=FALSE, verbose=TRUE, plotaic=TRUE, error_skip=FALSE, sample_shifts=FALSE, sample_threshold=2){
	
	n<-otree@nterm;nt<-dim(odata)[2]
	oldalphas<-sapply(oldfit,function(x)summary(x)$alpha)
	oldsigmas<-sapply(oldfit,function(x)summary(x)$sigma)
	oldparms<-data.frame(matrix(c(oldalphas,oldsigmas),nrow=2,byrow=T, dimnames=list(c("a","s"),names(odata))))
	odata2<-rbind(oldparms,odata)
	oldaic<-as.numeric(oldaic)
	
	uniqueshifts<-unique(oldshifts)
	k<-length(oldshifts);kk<-length(uniqueshifts)
	npair<-(kk*(kk-1))/2
	aics<-data.frame(X=rep(NA,npair),Y=rep(NA,npair),aic=rep(NA,npair),daic=rep(NA,npair))
	if(verbose)print(paste(kk, "unique regimes: testing", npair,"pairwise collapse(s)"),quote=F)
	counter<-0
for(i in 1:kk){
	for(j in 1:kk){
	if(i>j){
	counter<-counter+1
	tempshifts<-oldshifts
	tempshifts[tempshifts==uniqueshifts[j]]<-uniqueshifts[i]
	tempregs<-repaint(otree,regshifts=tempshifts)
	if(error_skip){
		te<-try( tempfit<-apply(odata2,2,function(x)hansen(x[-c(1,2)],otree,regimes=tempregs,sqrt.alpha=sqrt(x[1]),sigma=sqrt(x[2]))) )
		if(class(te)=="try-error"){
			print(paste("error fitting regimes",i,"and",j),quote=F)
			tempaic<-aic_threshold+1
		}else{
			tempLnL<-sum(sapply(tempfit,function(x)summary(x)$loglik))
			tempaic<-getAIC(tempLnL,k+nt*(2+kk-1),n*nt,TRUE)
		}
	}else{
	tempfit<-apply(odata2,2,function(x)hansen(x[-c(1,2)],otree,regimes=tempregs,sqrt.alpha=sqrt(x[1]),sigma=sqrt(x[2])))
	tempLnL<-sum(sapply(tempfit,function(x)summary(x)$loglik))
	tempaic<-getAIC(tempLnL,k+nt*(2+kk-1),n*nt,TRUE)
	}
	aics[counter,1:2]<-as.character(c(uniqueshifts[i],uniqueshifts[j]))
	aics[counter,3:4]<-c(tempaic,tempaic-oldaic)
	if(verbose)print(c(uniqueshifts[i],uniqueshifts[j],round(tempaic-oldaic,2)),quote=F)
}}}

newshifts<-oldshifts;best_aic<-oldaic;nreg<-kk;fit<-oldfit

if(any((aics$daic<(aic_threshold)))){

	df<-aics[which(aics$daic<(aic_threshold)),]
if(only_best==TRUE | dim(df)[1]==1){
	z<-which.min(df$daic)
	if(sample_shifts&df$daic[z]<(aic_threshold)){
		candidates<-which((df$daic-min(df$daic,na.rm=TRUE))<=sample_threshold&df$daic<(aic_threshold))
		if(verbose)print(paste("sampling 1 of",length(candidates), "models within",sample_threshold,"units of best AICc")) 
		if(length(candidates)>1)z<-sample(candidates,1)
	}
	if((dim(df)[1]==1)&verbose){print("1 cluster(s):",quote=F);print(2,quote=F)}
	if(verbose){print("collapsing:",quote=F);print(as.character(df[z,1:2]),quote=F)}
newshifts[which(newshifts==df$X[z])]<-rep(newshifts[which(newshifts==df$Y[z])][1],length(newshifts[which(newshifts==df$X[z])]))
	}else{

g<-igraph::graph.data.frame(df[,1:2],directed=FALSE)	
cl<-igraph::clusters(g,"weak")
vx<-unique(c(as.character(df$X),as.character(df$Y))) #same as V(g)
memb<-cl$memb
if(min(memb)==0)memb<-memb+1

if(verbose){print(paste(cl$no,"cluster(s):"),quote=F);print(cl$csize,quote=F)}

for(m in 1:cl$no){
	temp<-aics[aics$X%in%vx[memb==m] & aics$Y%in%vx[memb==m],]
if(all(temp$daic<aic_threshold,na.rm=TRUE)){ 
	v<-sort(vx[memb==m])
	}else{ 
	v<-sort(as.character(temp[which.min(temp$daic),1:2]))
	}
	if(verbose){print("collapsing:",quote=F);print(v,quote=F)}
	newshifts[which(newshifts%in%v[-1])]<-rep(v[1],length(newshifts[which(newshifts%in%v[-1])]))
	}
}
	regs<-repaint(otree,regshifts=newshifts)
	fit<-apply(odata2,2,function(x)hansen(x[-c(1,2)],otree,regimes=regs,sqrt.alpha=sqrt(x[1]),sigma=sqrt(x[2])))

	LnL<-sum(sapply(fit,function(x)summary(x)$loglik))
	nreg<-length(unique(newshifts))
	best_aic<-getAIC(LnL,k+nt*(2+nreg),n*nt,TRUE)
} 

xx<-summary(factor(newshifts));kk<-length(xx)
n_regimes<-c(k=k,kprime=kk,deltak=k-kk,c=sum(xx[xx>1]),kprime_conv=sum(xx>1),kprime_nonconv=sum(xx==1))

	if(plotaic){
		plot(aics$aic,ylim=c(min(best_aic,oldaic),oldaic+20),xlim=c(0,dim(aics)[1]));abline(h=oldaic-seq(0,30,by=2));abline(h=best_aic,col="red")
	}

return(list(fit=fit,all_aics=aics,aic=best_aic,savedshifts=newshifts,n_regimes=n_regimes))
}