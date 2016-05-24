addRegime<-
function(otree, odata, oldshifts, oldaic, oldfit, alloldaic=NULL, exclude=NULL, aic_threshold=0, verbose=FALSE, plotaic=FALSE, error_skip=FALSE, sample_shifts=FALSE, sample_threshold=2){
	if(is.null(oldshifts))oldshifts<-c("1"="a")
	Letters<-c(letters,paste("z",letters,sep=""),paste("zz",letters,sep="")) 
	Letters<-Letters[-which(Letters%in%oldshifts)]
	n<-otree@nterm;	nt<-dim(odata)[2]
	nodes<-otree@nodes
	uniqueshifts<-unique(oldshifts)
	k<-length(oldshifts)+1;kk<-length(uniqueshifts)+1
	aics<-LnLs<-rep(NA,length(nodes));names(aics)<-names(LnLs)<-nodes
	shifts<-character(length(nodes))
	fits<-list()
	if(class(oldfit)!="list")oldfit<-list(oldfit)
	
	oldalphas<-sapply(oldfit,function(x)summary(x)$alpha)
	oldsigmas<-sapply(oldfit,function(x)summary(x)$sigma)
	oldparms<-data.frame(matrix(c(oldalphas,oldsigmas),nrow=2,byrow=T, dimnames=list(c("a","s"),names(odata))))
	odata2<-rbind(oldparms,odata)

	if(!is.null(exclude)&!is.null(alloldaic)){
		old<-sort(alloldaic,decreasing=TRUE)
		Nexcluded<-floor(exclude*length(old))
		excluded<-names(old[0:Nexcluded])
			}else{		
		excluded<-NULL	
		}	
	skip<-c(excluded,names(oldshifts))

	if(verbose){
		print(paste("placing regime", k),quote=F)
		print(paste("testing ",length(nodes)-length(skip),"candidate models"),quote=F)
		}

for(i in 2:length(nodes)){
	if(i%in%skip==FALSE){
	shifts[i]<-Letters[1];names(shifts)[i]<-nodes[i]
	tempshifts<-c(oldshifts,Letters[1]);names(tempshifts)[k]<-i
	tempregs<-repaint(otree,regshifts=tempshifts)

	if(error_skip){
		te<-try( fits[[i]]<-apply(odata2,2,function(x)hansen(x[-c(1,2)],otree,regimes=tempregs,sqrt.alpha=sqrt(x[1]),sigma=sqrt(x[2]))) )
		if(class(te)=="try-error"){
			print(paste("error fitting regime",i),quote=F)
			LnLs[i]<-NA;aics[i]<-aic_threshold+9999
		}else{
			LnLs[i]<-sum(sapply(fits[[i]],function(x)summary(x)$loglik))
			aics[i]<-getAIC(LnLs[i],k+nt*(2+kk),n*nt,TRUE)
		}
	}else{
	fits[[i]]<-apply(odata2,2,function(x)hansen(x[-c(1,2)],otree,regimes=tempregs,sqrt.alpha=sqrt(x[1]),sigma=sqrt(x[2])))
	LnLs[i]<-sum(sapply(fits[[i]],function(x)summary(x)$loglik))
	aics[i]<-getAIC(LnLs[i],k+nt*(2+kk),n*nt,TRUE)
	}

	if(verbose)print(c(names(aics[i]),round(as.numeric(aics[i]-oldaic),2)), quote=F)
		}
	}
	best<-names(sort(aics))[1]
	if(sample_shifts&(aics[best]-oldaic)<(aic_threshold)){
		candidates<-aics[which((aics-min(aics,na.rm=TRUE))<=sample_threshold&(aics-oldaic)<(aic_threshold))]
		if(verbose)print(paste("sampling 1 of",length(candidates), "models within",sample_threshold,"units of best AICc")) 
		if(length(candidates)>1)best<-names(sample(candidates,1))
	}
	newshifts<-c(oldshifts,shifts[as.numeric(best)])
	xx<-summary(factor(newshifts))
	n_regimes<-c(k=k,kprime=kk,deltak=k-kk,c=sum(xx[xx>1]),kprime_conv=sum(xx>1),kprime_nonconv=sum(xx==1))
	if(plotaic){
		plot(aics,ylim=c(min(oldaic,aics[best]),oldaic+20),xlim=c(0,dim(odata)[1]),main=k);abline(h=oldaic-seq(0,30,by=2));abline(h=aics[best],col="red")
	}
	if(verbose){
		print(paste("old AIC =",round(oldaic,2)),quote=F)
		print(paste("new AIC =",round(as.numeric(aics[best]),2)),quote=F)
		if((aics[best]-oldaic)<(aic_threshold))print(paste("adding regime shift at node",names(aics[best])),quote=F)
		}
return(list(fit=fits[[as.numeric(best)]],all_aic=aics,aic=aics[best],savedshifts=newshifts,n_regimes=n_regimes))	}
