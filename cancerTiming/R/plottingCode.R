
plotAlleleByPosition<-function(mutData,segmentData=NULL,whChr=1:22,chromosomeId="chr",sampleName=NULL,sample=FALSE,
	tumorAFId,positionId,type="mutation",
    startId="start",endId="end",segFactorId,tCNId,MarkId,segColors,normCont=NULL,
	addData=NULL,addColor="red",col="black",pch=1,lwd=2,xlim,ylab="Allele Frequency",...){
	
	givenNames<- c(chromosomeId,tumorAFId,positionId)
	changeToNames<-c("chr","tAF","pos")
	mToNames<-match(givenNames,names(mutData))
	if(any(is.na(mToNames))) stop(givenNames[is.na(mToNames)],"not found in names of mutData")
	names(mutData)[mToNames]<-changeToNames
	mutData$chr<-as.character(mutData$chr)
	whChr<-as.character(whChr)
	doLOH<-!is.null(addData)
	doSeg<-!is.null(segmentData)
    doExpAF<-!is.null(normCont) & doSeg
    doMarkLOHSeg<-!missing(MarkId)
	if(doLOH){
		mToNames<-match(givenNames,names(addData))
		if(any(is.na(mToNames))) stop(givenNames[is.na(mToNames)],"not found in names of addData")
		names(addData)[mToNames]<-changeToNames
		addData$chr<-as.character(addData$chr)
		
	}
	if(doSeg){
		givenNames<-c(chromosomeId,startId,endId,segFactorId)
        changeToNames<-c(changeToNames[1],"start","end","segFactor")
        if(doExpAF){
            if(tCNId ==segFactorId){
               segmentData$totalCopyNumber<-segmentData[[segFactorId]]
            }
            else{
                givenNames<-c(givenNames,tCNId)
                changeToNames<-c(changeToNames,"totalCopyNumber")
            }
         }
        if(doMarkLOHSeg){
                  givenNames<-c(givenNames,MarkId)
                changeToNames<-c(changeToNames,"LOH")
     
        }
		mToNames<-match(givenNames,names(segmentData))
		if(any(is.na(mToNames))) stop(givenNames[is.na(mToNames)],"not found in names of segmentData")
		names(segmentData)[mToNames]<-changeToNames
		segmentData$chr<-as.character(segmentData$chr)
		if(!is.factor(segmentData$segFactor)) segmentData$segFactor<-factor(segmentData$segFactor)
		if(missing(segColors)) segColors<-rep(palette()[-1],length=nlevels(segmentData$segFactor))
		if(length(segColors)!=nlevels(segmentData$segFactor)) stop("Invalid length of segColors -- doesn't match number of levels of segFactor")
        names(segColors)<-levels(segmentData$segFactor)
		
	}
	if(!any(whChr %in% mutData$chr)) warning("No values of whChr found in mutData")
	for( i in whChr){
	    chrdataSNP <- mutData[mutData$chr==i,]#subset(mutData, chr == i);
        if(sample){
            if(nrow(chrdataSNP)>10000){
                ind<-sample(1:nrow(chrdataSNP),size=10000,replace=FALSE)
                chrdataSNP<-chrdataSNP[ind,]
                }
        }
		if(missing(xlim)){
			xlim<-range(chrdataSNP$pos)
			if(doSeg) xlim<-range(c(xlim,segmentData$start[segmentData$chr==i],segmentData$end[segmentData$chr==i]))
		}
	    if(doLOH) chrdataLOH <- addData[addData$chr==i,]#subset(addData, chr == i );
	    if(nrow(chrdataSNP)==0) next;
		if(doLOH && nrow(chrdataLOH) > 0){
			with(chrdataLOH, plot(pos, tAF, ylim=c(0, 1), col=addColor,xlim=xlim,ylab=ylab,...))
	    	with(chrdataSNP, points(pos, tAF,col=col,pch=pch))
		}
		else{
			with(chrdataSNP, plot(pos, tAF, ylim=c(0, 1), xlim=xlim,ylab=ylab,col=col,pch=pch,...))
		}
		axis(3,at=par("usr")[1],paste(sampleName,": Chromosome", i))
		if(doSeg){ #put rectangle on those think are cnloh
			if(i %in% segmentData$chr){
				whchr<-which(segmentData$chr==i)
				for(ind in whchr){
					col<-segColors[as.character(segmentData$segFactor[[ind]])]
					pct<-0.025*(diff(par("usr")[3:4]))
					
                    dens<-if(doMarkLOHSeg && segmentData$LOH[[ind]]) 10 else NULL
					rect(xleft=segmentData$start[ind], ybottom=par("usr")[3]-2*pct, xright=segmentData$end[ind],lwd=2, ytop=par("usr")[3], col = col,xpd=NA, density=dens)
                    if(doExpAF){
                        if(!is.na(segmentData$totalCopyNumber[ind])){
                            if(type=="mutation"){
                                AF<-allAF(segmentData$totalCopyNumber[ind],normCont,type="mutation")
                                if(length(AF[[1]])>0) segments(x0=segmentData$start[ind],x1=segmentData$end[ind],y0=AF[[1]],y1=AF[[1]],col=col,lwd=lwd)
                            }
                            else{
                                 AFHet<-allAF(segmentData$totalCopyNumber[ind],normCont,type="SNPHet")
                                 AFHomo<-allAF(segmentData$totalCopyNumber[ind],normCont,type="SNPHomo")
                                if(length(AFHet[[1]])>0) segments(x0=segmentData$start[ind],x1=segmentData$end[ind],y0=AFHet[[1]],y1=AFHet[[1]],col=col,lwd=lwd)
                                if(length(AFHomo[[1]])>0) segments(x0=segmentData$start[ind],x1=segmentData$end[ind],y0=AFHomo[[1]],y1=AFHomo[[1]],col=col,lty=3,lwd=lwd)
                           
                            }
                        }
                    }
				}
			}
		}
	}
	invisible(segColors)
}

#for total number of copies of tumor and proportion of normal contamination, calculated possible expected allele frequencies
plotAlleleDensity<-function(af,depth,groupingId,totalCopy,groupCol=palette(),normCont=0,type="mutation",minDepth=40,lineCols=c("grey","tan4"),minMut=40,histogram=FALSE){
	if(missing(groupingId)) groupingId<-rep("All Entries",length=length(af))
	ids<-as.character(groupingId)
	#for id, make density/histogram of the allele frequencies (af)
	doLines<-!missing(normCont)
	if(doLines) lineCols<-rep(lineCols,length=length(normCont))

	indsList<-tapply(1:length(ids),ids,function(x){x})
	x<-by(cbind(af=af,depth=depth),list(ids),function(x){x$af[x$depth>minDepth]},simplify=FALSE)
	totdepth<-by(cbind(af=af,depth=depth),list(ids),function(x){x$depth[x$depth>minDepth]},simplify=FALSE)
	col<-rep(groupCol,length=length(x))	
	whPlot<-sapply(x,length)>minMut
	if(sum(whPlot)>0){
		x<-x[whPlot]
		col<-col[whPlot]
		totdepth<-totdepth[whPlot]
		#ids<-ids[whPlot]
		#calculate the percent tumor for this total number of copies (assume normal diploid):
		if(is.null(names(normCont))){
			names(normCont)<-paste("Norm.Cont. Est",1:length(normCont))
		}
		pTumor<-sapply(normCont,function(x){totalCopy*(1-x)/(totalCopy*(1-x)+2*(x))})
		lineFun<-function(total){
                AF<-allAF(total,normCont,type=type)
                ltys<-rep(1:3,length=total)
                mapply(AF,lineCols,FUN=function(x,col){
					abline(v=x,lty=ltys,col=col,lwd=2)					
				})
				legend("topleft",c(names(AF[[1]]),paste(names(normCont),"(",names(AF),")",sep="")),
					lty=c(ltys,rep(1,length(AF))),col=c(rep("black",length=total),lineCols),lwd=2)
				
		}
		if(histogram){#plot individual histograms
			mapply(x,totdepth,names(x),FUN=function(dat,depth,id){
				hist(dat,xlim=c(0,1),xlab="Observed Allele Frequency",main=paste("Region",id),
				sub=paste(length(dat),"mutations with minimum depth",min(depth),"and max depth",max(depth))
				)
			if(doLines) lineFun(totalCopy)
			})
			
			
		}
		else{
			multidensity(x,col=col,lwd=1,lty=1,xlim=c(0,1),ylab="Density",xlab="Observed Allele Frequency",
				main=paste("Total Copy number=",totalCopy),
				sub=paste("Only regions with at least",minMut,"mutations that each have read depth of at least",minDepth))
			if(doLines) lineFun(totalCopy)
		}
	}
	invisible(x)
}
multidensity<-function(x,col=palette(),lwd=1,lty=1,xlim,ylab="Density",...){
	if(missing(xlim)){
		xvals<-unlist(lapply(x,function(z){range(z[is.finite(z)],na.rm=TRUE)}))
		xlim<-range(xvals[is.finite(xvals)])
	}
	dens<-lapply(x,function(x){
		density(x[is.finite(x)])
	})
	yvals<-unlist(lapply(dens,function(x){x$y}))
	ylim<-range(yvals[is.finite(yvals)])
	plot(0,type="n",xlim=xlim,ylim=ylim,ylab=ylab,...)
	out<-mapply(dens,rep(col,length=length(dens)),rep(lwd,length=length(dens)),rep(lty,length=length(dens)),FUN=function(x,col,lwd,lty){lines(x,col=col,lwd=lwd,lty=lty)})
}

plotCopies<-function(x,y,nX,nY,xleg,yleg,onlyPositive=TRUE,equalAxis=TRUE,integerLegend=TRUE,xlim,ylim,...){
	if(onlyPositive){
		whneg<- which(x<0 | y<0)
		if(length(whneg)>0){
			x<-x[-whneg]
			y<-y[-whneg]
			nX<-nX[-whneg]
			nY<-nY[-whneg]
		}
		if(length(x)==0) stop("no positive values in x or y")		
	}
	makeValues<-function(x,parValues){
		names(x)<-x
		if(onlyPositive){
			x[x<0]<-NA
			noMatch<-any(is.na(x))
			if(!noMatch) parValues<-parValues[-1]
			names(x)[is.na(x)]<-"No clear correspondence"
		}
		else parValues<-parValues[-1]
		numValues<-sort(unique(x[!is.na(x)]))
		maxVal<-max(x,na.rm=TRUE)
		xFac<-if(integerLegend) factor(names(x),levels=if(noMatch) c("No clear correspondence",0:maxVal) else 0:maxVal) else factor(if(is.numeric(x)) round(x,3) else x)
		parValuesVec<-rep(parValues,nlevels(xFac))[xFac]
		return(list(parValuesVec=parValuesVec,parValues=parValues,factor=xFac))
	}
	pchC<-makeValues(nX,c(1,c(22,24,23,25)))
	colC<-makeValues(nY,palette()[-6])		
	if(equalAxis){
		lim<-if(onlyPositive) range(c(0,x,y),na.rm=TRUE) else range(c(x,y),na.rm=TRUE)
		xlim<-ylim<-lim
	}
	else{
		if(missing(xlim)) xlim<-range(x,na.rm=TRUE)
		if(missing(ylim)) ylim<-range(y,na.rm=TRUE)
	}
	plot(x,y,col="black",xlim=xlim,ylim=ylim,bg=colC$parValuesVec,pch=pchC$parValuesVec,...)
	if(missing(yleg)) yleg<-""
	if(missing(xleg)) xleg<-""
	legend("topleft",levels(colC$factor),fill=colC$parValues,title=yleg,bty="n")
	legend("bottomright",levels(pchC$factor),pch=pchC$parValues,pt.bg="black",title=xleg,bty="n")
	invisible(list(col=colC,pch=pchC))
}

plotSegmentation<-function(segs,valId,col=palette(),lty=1,lwd=2,xlim,ylim,xlab="Position",ylab=valId,...){
	if(missing(col)) col<-rep(col,length=length(segs))
	col<-rep(col,length=length(segs))
	lty<-rep(lty,length=length(segs))
	names(col)<-names(segs)
	if(missing(xlim)){
		x<-unlist(lapply(segs,function(x){c(x$start,x$end)}))
		xlim<-range(x,na.rm=TRUE)
	}
	if(missing(ylim)){
		y<-unlist(lapply(segs,function(x){unlist(x[,valId])}))
		ylim<-range(y,na.rm=TRUE)
	} 
	plot(0,type="n",xlim=xlim,ylim=ylim,xlab=xlab,ylab=ylab,...)
	for(ii in 1:length(segs)){
		seg<-segs[[ii]]
		if(nrow(seg)>0) segments(x0=seg$start,x1=seg$end,y0=seg[[valId]],lty=lty[[ii]],col=col[[ii]],lwd=lwd)			
	}
	invisible(list(col=col,lty=lty))
}

#segs should be factor saying what segment in
plotSeqCount<-function(position,t_count,n_count,ylim=NULL,normFac=1,segs,segColors=palette(),...){
	y<-log2(t_count/n_count)-log2(normFac)
	wh<-which(!is.na(y) & is.finite(y))
	x<-position[wh]
	y<-y[wh]
	if(is.null(ylim)) ylim<-quantile(y[wh],probs=c(.05,.95),na.rm=TRUE)
	plot(x,y,ylim=ylim,...)
	#lines(lowess(x=x,y=y), col = "red")
	if(!missing(segs)){
		segs<-segs[wh]
		segAve<-tapply(y,segs,mean)
		segBeg<-tapply(x,segs,min)
		segEnd<-tapply(x,segs,max)	
		segments(x0=segBeg,x1=segEnd,y0=segAve,col=segColors,lwd=2)
	}
	abline(h=0)
	invisible(list(x=x,y=y))
	
}