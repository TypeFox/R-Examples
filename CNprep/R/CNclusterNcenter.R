CNclusterNcenter <-
function(segrat,blsize,minjoin,ntrial,bestbic,
	modelNames,cweight,bstimes,chromrange,seedme){
	.lec.CreateStream(segrat$stream)
        .lec.SetSeed(segrat$stream,seedme)
        for( j in 1:segrat$sub).lec.ResetNextSubstream(segrat$stream)
        .lec.CurrentStream(segrat$stream)
	startcol<-"StartProbe"
	endcol<-"EndProbe"
	chromcol<-"chrom"
	medcol<-"segmedian"
	madcol<-"segmad"
	segrat$seg<-cbind(segrat$seg,t(apply(segrat$seg[,c(startcol,endcol,chromcol),
		drop=F],1,smedmad,v=segrat$rat)))
	dimnames(segrat$seg)[[2]]<-c(startcol,endcol,chromcol,medcol,madcol)
	seguse<-segrat$seg[segrat$seg[,chromcol]%in%chromrange,,drop=F]
	aux<-rep(0,length(segrat$rat))
	aux[seguse[,startcol]]<-1
	aux[seguse[,endcol]]<-(-1)
	aux<-cumsum(aux)
	aux[seguse[,endcol]]<-1
	ratuse<-segrat$rat[aux==1]
	for(j in 1:ntrial){
		aaa<-segsample(seguse,ratuse,blocksize=blsize)
		if(all(unique(aaa[,3])==0)){ aaa[,3]<-1e-10 }
		emfit<-Mclust(aaa[,3],maxG=15,modelNames=modelNames)
		if(emfit$bic>=bestbic){
			bestaaa<-aaa
			bestem<-emfit
			bestbic<-emfit$bic
		}
	}
	newem<-consolidate(bestem,minjoin)
	newem<-get.center(newem,cweight)
	if(length(bestem$parameters$mean)==1){ profcenter<-median(bestaaa[,3])
        }else{ profcenter<-weighted.median(bestaaa[,3],newem$z[,newem$center]) }
	mediandev<-segrat$seg[,medcol]-profcenter
	segs<-segsample(segrat$seg,segrat$rat,times=bstimes)
	if(all(unique(aaa[,3])==1e-10)){ segs[segs[,3]==0,3]<-1e-10 }
	segzall<-getz(segs[,3],bestem,newem$groups,times=bstimes)
	centerz<-segzall[,newem$center]
	maxz<-segzall[nrow(segzall)*(max.col(segzall)-1)+(1:nrow(segzall))]
	maxzcol<-max.col(segzall)
	maxzmean<-newem$mu[maxzcol]-newem$mu[newem$center]
	maxzsigma<-sqrt(newem$sigmasq[maxzcol])
	cpb<-centerprob(segs[,3],bestem,newem$groups,times=bstimes,newem$center)
	w<-t(matrix(nrow=bstimes,data=segs[,3]))
	segerr<-sqrt(apply(w,1,var,na.rm=T))
	.lec.CurrentStreamEnd()
        .lec.DeleteStream(segrat$stream)
	return(cbind(segrat$seg[,medcol],segrat$seg[,madcol],mediandev,segerr,centerz,
		cpb,maxz,maxzmean,maxzsigma))
}
