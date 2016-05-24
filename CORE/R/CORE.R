CORE <-
function(dataIn,keep=NULL,startcol="start",endcol="end",
chromcol="chrom",weightcol="weight",maxmark=1,minscore=0,pow=1,
assoc=c("I","J","P"),nshuffle=0,boundaries=NULL,seedme=sample(1e8,1),
shufflemethod=c("SIMPLE","RESCALE"),tiny=-1,
distrib=c("vanilla","Rparallel","Grid"),njobs=1,qmem=NA){
	doshuffles<-"NO"
	shufflemethod<-match.arg(shufflemethod)
	assoc<-match.arg(assoc)
	if(!(class(dataIn)%in%c("CORE","matrix","data.frame")))
		stop("invalid class of argument dataIn")
	if(class(dataIn)=="CORE")if(((!("assoc"%in%keep))&(dataIn$assoc!=assoc))|
		((!("pow"%in%keep))&(dataIn$pow!=pow))){
		cankeep<-c("maxmark","minscore","nshuffle","boundaries","seedme",
				"shufflemethod","tiny")
		for(item in intersect(keep,cankeep))assign(item,dataIn[[item]])
		dataIn<-dataIn$input
		startcol<-"start"
		endcol<-"end"
		chromcol<-"chrom"
		weightcol<-"weight"
	}
	if(class(dataIn)!="CORE"){	#get cores from scratch
		if(!(startcol%in%dimnames(dataIn)[[2]]))
			stop("start column missing in input data")
		if(!(endcol%in%dimnames(dataIn)[[2]]))
			stop("end column missing in input data")
		if(!(chromcol%in%dimnames(dataIn)[[2]])){
			dataIn<-cbind(dataIn,rep(1,nrow(dataIn)))
			dimnames(dataIn)[[2]][ncol(dataIn)]<-chromcol
		}
		if(!(weightcol%in%dimnames(dataIn)[[2]])){
			dataIn<-cbind(dataIn,rep(1,nrow(dataIn)))
			dimnames(dataIn)[[2]][ncol(dataIn)]<-weightcol
		}
		z<-as.matrix(dataIn[,c(chromcol,startcol,endcol,weightcol),drop=F])
		#rm(dataIn)
		dimnames(z)[[2]]<-c("chrom","start","end","weight")
		if(!is.null(boundaries)){
			if(!(startcol%in%dimnames(boundaries)[[2]]))
				stop("start column missing in boundary table")
			if(!(endcol%in%dimnames(boundaries)[[2]]))
				stop("end column missing in boundary table")
			if(!(chromcol%in%dimnames(boundaries)[[2]]))
				stop("chrom column missing in boundary table")
			boundaries<-as.matrix(boundaries[,c(chromcol,startcol,endcol),drop=F])
		}
		result<-switch(assoc,
			I=ICOREiteration(z,maxmark,pow,minscore),
			J=JCOREiteration(z,maxmark,pow,minscore),
			P=PCOREiteration(z,maxmark,minscore)
		)
		returnme<-list(input=z,call=match.call(),minscore=minscore,maxmark=maxmark,
			pow=pow,assoc=assoc,coreTable=result,seedme=seedme,boundaries=boundaries,
			shufflemethod=shufflemethod,nshuffle=nshuffle,tiny=tiny)
		if(nshuffle>0)doshuffles<-"FROMSCRATCH"
	}
	#If dataIn is a CORE, determine whether computation is to be continued
	else{ 
		assoc<-dataIn$assoc
		pow<-dataIn$pow
		cankeep<-c("maxmark","minscore",
		"nshuffle","boundaries","seedme","shufflemethod","tiny")
		for(item in intersect(keep,cankeep))assign(item,dataIn[[item]])
		returnme<-list(input=dataIn$input,call=match.call(),minscore=minscore,
			maxmark=maxmark,pow=dataIn$pow,assoc=dataIn$assoc,
			shufflemethod=shufflemethod,seedme=seedme,nshuffle=nshuffle,tiny=tiny,
			boundaries=boundaries)
		if(dataIn$coreTable[nrow(dataIn$coreTable),"score"]>=minscore&
			nrow(dataIn$coreTable)<maxmark){
			z<-dataIn$input
			if(assoc=="I"){
				for(i in 1:nrow(dataIn$coreTable)){
					fixus<-z[,"chrom"]==dataIn$coreTable[i,"chrom"]&
						z[,"start"]<=dataIn$coreTable[i,"start"]&
						z[,"end"]>=dataIn$coreTable[i,"end"]
					z[fixus,"weight"]<-z[fixus,"weight"]*
						(1-((dataIn$coreTable[i,"end"]-dataIn$coreTable[i,"start"]+1)/
						(z[fixus,"end"]-z[fixus,"start"]+1))^dataIn$pow)
				}
				z<-z[z[,"weight"]>=tiny,,drop=F]
				result<-ICOREiteration(z,maxmark-nrow(dataIn$coreTable),
					dataIn$pow,minscore)
			}
			if(assoc=="J"){
				for(i in 1:nrow(dataIn$coreTable)){
					fixus<-z[,"chrom"]==dataIn$coreTable[i,"chrom"]&
						pmax(z[,"start"],dataIn$coreTable[i,"start"])<=
						pmin(z[,"end"],dataIn$coreTable[i,"end"])
					z[fixus,"weight"]<-z[fixus,"weight"]*
						(1-((pmin(dataIn$coreTable[i,"end"],z[fixus,"end"])-
						pmax(dataIn$coreTable[i,"start"],z[fixus,"start"])+1)/
						(pmax(dataIn$coreTable[i,"end"],z[fixus,"end"])-
						pmin(dataIn$coreTable[i,"start"],z[fixus,"start"])+1))^
						dataIn$pow)
				}
				z<-z[z[,"weight"]>=tiny,,drop=F]
				result<-JCOREiteration(z,maxmark-nrow(dataIn$coreTable),
					dataIn$pow,minscore)
			}
			if(assoc=="P"){
				for(i in 1:nrow(dataIn$coreTable)){
					fixus<-z[,"chrom"]==dataIn$coreTable[i,"chrom"]&
						z[,"start"]<=dataIn$coreTable[i,"start"]&
						z[,"end"]>=dataIn$coreTable[i,"end"]
					z<-z[-fixus,,drop=F]
				}
				result<-PCOREiteration(z,maxmark-nrow(dataIn$coreTable),
					minscore)
			}
			returnme$coreTable<-rbind(dataIn$coreTable,result)
		}
		else returnme$coreTable<-dataIn$coreTable[which(dataIn$coreTable[1:
			min(maxmark,nrow(dataIn$coreTable)),"score"]>=minscore),,drop=F]
		if(nshuffle>0){
			if(dataIn$shufflemethod!=shufflemethod|dataIn$seedme!=seedme|
				nrow(dataIn$coreTable)<nrow(returnme$coreTable)|dataIn$nshuffle==0)
				doshuffles<-"FROMSCRATCH"
			else if(dataIn$nshuffle<nshuffle)doshuffles<-"ADD"
			else returnme$simscores<-
				dataIn$simscores[1:nrow(returnme$coreTable),1:nshuffle,drop=F]
		}
	}
	class(returnme)<-"CORE"
	if(is.null(boundaries)){
		y<-returnme$input[order(returnme$input[,"chrom"]),,drop=F]
		boundaries<-cbind(unique(y[,"chrom"]),
			tapply(X=y[,"start"],INDEX=y[,"chrom"],FUN=min),
			tapply(X=y[,"end"],INDEX=y[,"chrom"],FUN=max))
	}
	dimnames(boundaries)[[2]]<-c("chrom","start","end")
	if(assoc=="I") randfun<-ICORErandomized
	if(assoc=="J") randfun<-JCORErandomized
	if(assoc=="P") randfun<-PCORErandomized
	distrib<-match.arg(distrib)
	returnme<-Rparallel(randfun,distrib,doshuffles,nshuffle,dataIn,returnme,boundaries,njobs,qmem)
	return(returnme)
}
