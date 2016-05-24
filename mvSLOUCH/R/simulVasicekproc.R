.ouch.simulate<-function(step,duration,modelParams,regimes=NULL,regimes.times=NULL,mCov=NULL){
    if (is.null(regimes)){
	if (is.null(colnames(modelParams$mPsi))){
	    regimes<-c("1")
	    colnames(modelParams$mPsi)<-as.character(1:ncol(modelParams$mPsi))
	}else{regimes<-colnames(modelParams$mPsi)[1]}
    }
    if (is.null(regimes.times)){regimes.times<-seq(0,duration,length.out=length(regimes)+1)}

    kY<-nrow(modelParams$A)
    modelParams$precalcMatrices<-.decompEigenA.S(modelParams,NULL,NULL,list(bCalcA=TRUE,bCovCalc=TRUE,dzetacalc=FALSE,lexptcalc=FALSE,kappacalc=FALSE,interceptcalc=FALSE),NULL)
    mTrajectory<-matrix(NA,ncol=kY+1,nrow=ceiling(duration/step)+1)
    mTrajectory[1,]<-c(0,c(modelParams$vY0[,1]))
    i<-2
    curr.time<-step
    lRegs<-list(times=NULL,exptjA=NULL,regimes=NULL)
    while(curr.time<duration){
	if (!is.null(regimes)){
	    vRegTimesIndex<-intersect(which(regimes.times>curr.time-step),which(regimes.times<=curr.time))
	    vRegTimesIndex<-c(max(which(regimes.times<=curr.time-step)),vRegTimesIndex)
	    lRegs$times<-regimes.times[vRegTimesIndex]
	    lRegs$times[1]<-curr.time-step
	    if (lRegs$times[length(lRegs$times)]<curr.time){lRegs$times<-c(lRegs$times,curr.time)}
	    lRegs$regimes<-rep(NA,length(lRegs$times)-1)
	    for (j in 1:length(lRegs$regimes)){
		if (lRegs$times[j+1]>=max(regimes.times)){lRegs$regimes[j]<-regimes[length(regimes.times)]}
		else{lRegs$regimes[j]<-regimes[min(which(regimes.times>=lRegs$times[j+1]))-1]}
	    }
	    lRegs$times<-lRegs$times-(curr.time-step)
	}
	mTrajectory[i,]<-c(curr.time,.draw.ouch(step,mTrajectory[i-1,2:(kY+1)],modelParams,mCov,lRegs))
	i<-i+1
	curr.time<-curr.time+step
    }
    mTrajectory<-Re(mTrajectory)
    mTrajectory
}


.draw.ouch<-function(curr.time,Y,modelParams,mCov=NULL,lRegs=list(times=NULL,exptjA=NULL,regimes=NULL)){
    expmtA<-modelParams$precalcMatrices[[1]]$eigA$vectors%*%diag(exp(-modelParams$precalcMatrices[[1]]$eigA$values*curr.time),length(Y),length(Y))%*%modelParams$precalcMatrices[[1]]$invP
    if (!is.null(lRegs$times)){lRegs$exptjA<-sapply(lRegs$times,function(t1,eigA,invP){eigA$vectors%*%diag(exp(eigA$values*t1),length(eigA$values),length(eigA$values))%*%invP},eigA=modelParams$precalcMatrices[[1]]$eigA,invP=modelParams$precalcMatrices[[1]]$invP,simplify=FALSE)}
    modelParams$vY0<-Y
    vMean<-.calc.mean.ouch.mv(expmtA,Y,modelParams$mPsi,modelParams$mPsi0,lRegs$exptjA,lRegs$regimes)
    if (is.null(mCov)){mCov<-.calc.cov.ouch.mv(curr.time,modelParams$precalcMatrices[[1]],modelParams$precalcMatrices[[2]])}
    rmvnorm(n=1,mean=vMean,sigma=mCov)
}

.bm.simulate<-function(step,duration,modelParams,regimes=NULL,regimes.times=NULL,mCov=NULL){
    SigmaSq<-modelParams$Sxx%*%t(modelParams$Sxx)
    k<-nrow(SigmaSq)
    vX0<-modelParams$vX0
    
    mTrajectory<-matrix(0,ncol=k+1,nrow=ceiling(duration/step)+1)
    mTrajectory[1,]<-c(0,vX0)
    i<-2
    curr.time<-step
    while(curr.time<duration){
	mTrajectory[i,]<-c(curr.time,.draw.bm(step,mTrajectory[i-1,2:(k+1)],SigmaSq))
	i<-i+1
	curr.time<-curr.time+step
    }
    mTrajectory
}

.draw.bm<-function(curr.time,Y,Sigmasq){
    rmvnorm(n=1,mean=Y,sigma=curr.time*Sigmasq)
}

.mvslouch.simulate<-function(step,duration,modelParams,regimes=NULL,regimes.times=NULL,mCov=NULL){
    if (is.null(regimes)){
	if (is.null(colnames(modelParams$mPsi))){
	    regimes<-c("1")
	    colnames(modelParams$mPsi)<-as.character(1:ncol(modelParams$mPsi))
	}else{regimes<-colnames(modelParams$mPsi)[1]}
    }
    if (is.null(regimes.times)){regimes.times<-seq(0,1,length.out=length(regimes)+1)}

    kY<-nrow(modelParams$A)
    kX<-ncol(modelParams$B)
    modelParams$precalcMatrices<-.decompEigenA.S(modelParams,NULL,NULL,list(bCalcA=TRUE,bCovCalc=TRUE,dzetacalc=FALSE,lexptcalc=FALSE,kappacalc=FALSE,interceptcalc=FALSE),NULL)
    mTrajectory<-matrix(NA,ncol=kY+kX+1,nrow=ceiling(duration/step)+1)
    mTrajectory[1,]<-c(0,c(modelParams$vY0[,1],modelParams$vX0[,1]))
    i<-2
    curr.time<-step
    lRegs<-list(times=NULL,exptjA=NULL,regimes=NULL)
    while(curr.time<duration){
	if (!is.null(regimes)){
	    vRegTimesIndex<-intersect(which(regimes.times>curr.time-step),which(regimes.times<=curr.time))
	    vRegTimesIndex<-c(max(which(regimes.times<=curr.time-step)),vRegTimesIndex)
	    lRegs$times<-regimes.times[vRegTimesIndex]
	    lRegs$times[1]<-curr.time-step
	    if (lRegs$times[length(lRegs$times)]<curr.time){lRegs$times<-c(lRegs$times,curr.time)}
	    lRegs$regimes<-rep(NA,length(lRegs$times)-1)
	    for (j in 1:length(lRegs$regimes)){

		if (lRegs$times[j+1]>=max(regimes.times)){lRegs$regimes[j]<-regimes[length(regimes.times)]}
		else{lRegs$regimes[j]<-regimes[min(which(regimes.times>=lRegs$times[j+1]))-1]}
	    }
	    lRegs$times<-lRegs$times-(curr.time-step)
	}
	mTrajectory[i,]<-c(curr.time,.draw.mvslouch(step,mTrajectory[i-1,2:(kY+1)],mTrajectory[i-1,(kY+2):(kY+kX+1)],modelParams,mCov,lRegs))
	i<-i+1
	curr.time<-curr.time+step
    }
    mTrajectory<-Re(mTrajectory)
    mTrajectory
}

.draw.mvslouch<-function(curr.time,Y,vX,modelParams,mCov=NULL,lRegs=list(times=NULL,exptjA=NULL,regimes=NULL)){
    expmtA<-modelParams$precalcMatrices[[1]]$eigA$vectors%*%diag(exp(-modelParams$precalcMatrices[[1]]$eigA$values*curr.time),length(Y),length(Y))%*%modelParams$precalcMatrices[[1]]$invP
    if (!is.null(lRegs$times)){lRegs$exptjA<-sapply(lRegs$times,function(t1,eigA,invP){eigA$vectors%*%diag(exp(eigA$values*t1),length(eigA$values),length(eigA$values))%*%invP},eigA=modelParams$precalcMatrices[[1]]$eigA,invP=modelParams$precalcMatrices[[1]]$invP,simplify=FALSE)}
    modelParams$vY0<-Y
    modelParams$vX0<-vX
    vMean<-.calc.mean.slouch.mv(expmtA,modelParams$precalcMatrices[[1]]$A1B,Y,vX,modelParams$mPsi,modelParams$mPsi0,lRegs$exptjA,lRegs$regimes)
    if (is.null(mCov)){mCov<-.calc.cov.slouch.mv(curr.time,modelParams$precalcMatrices[[1]],modelParams$precalcMatrices[[2]])}
    rmvnorm(n=1,mean=vMean,sigma=mCov)
}


.drawOneLineage<-function(mTraject,vColours,kY,kX,sText="",newplot=TRUE,boptDraw=FALSE,legendplace="bottomleft"){
    if(boptDraw){vcolsplot<-2:(2*kY+kX+1)}else{vcolsplot<-2:(kY+kX+1)}
    if (newplot&&is.null(colnames(mTraject))){
	colnames(mTraject)<-1:ncol(mTraject)
    	colnames(mTraject)[1]<-"time"
    	for(i in vcolsplot){
    	    if (i<=kY+1){colnames(mTraject)[i]<-paste("Y",i-1)}
	    else{
		if(i<=kY+kX+1){colnames(mTraject)[i]<-paste("T",i-kY-1)}
		else{colnames(mTraject)[i]<-paste("optY",i-kY-kX-1)}
	    }
	}
    }
    if(newplot){
        plot(rep(mTraject[,1],length(vcolsplot)),c(mTraject[,vcolsplot]),col="white",pch=19,cex=0.3,main=sText,xlab=colnames(mTraject)[1],ylab="(Y,X)")
	legendText<-rep("",length(vcolsplot))
    }    
    for(i in vcolsplot){
	points(mTraject[,1],mTraject[,i],pch=19,cex=0.3,col=vColours[i-1])
	if (newplot){legendText[i-1]<-colnames(mTraject)[i]}	
    }
    if(newplot&&!is.na(legendplace)){legend(legendplace,legend=legendText,col=vColours,pch=19)}    
}

