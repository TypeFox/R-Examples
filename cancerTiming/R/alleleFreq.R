# #code to skeleton for documentation; DO NOT RUN AGAIN, will overwrite the changes!
# source("~/Documents/RfunctionsGenerally/InternalRPackages/cancerTiming/R/eventTiming.R", chdir = FALSE)
# source("~/Documents/RfunctionsGenerally/InternalRPackages/cancerTiming/R/alleleFreq.R", chdir = FALSE)
# source("~/Documents/RfunctionsGenerally/InternalRPackages/cancerTiming/R/bootstrapEventTiming.R", chdir = FALSE)
# sapply(ls(),function(x){prompt(object=get(x),name=x,filename=paste("~/Documents/RfunctionsGenerally/InternalRPackages/cancerTiming/man/",x,".Rd",sep=""))})

.checkLength<-function(x,nam,input,inputName){if(length(x)>1){if(length(x)!=length(input)) stop(paste(nam,"must be of same length as",inputName)); return(x)} else return(rep(x,length(input)))}
contAF<-function(trueAF,totalCopy,normCont=0,totalCopyNormal=2,type=c("mutation","SNPHet","SNPHomo")){
	type<-match.arg(type)
	if(any(normCont>1 | normCont<0)) stop("normCont must be between 0 and 1")
	totalCopy<-.checkLength(totalCopy,"totalCopy",trueAF,"trueAF")
	normCont<-.checkLength(normCont,"normCont",trueAF,"trueAF")
	totalCopyNormal<-.checkLength(totalCopyNormal,"totalCopyNormal",trueAF,"trueAF")
	const<-switch(type,"mutation"=0,"SNPHet"=1,SNPHomo=totalCopyNormal) #number of copies in the normal that contain the variant
	af<-(trueAF*totalCopy*(1-normCont)+const*normCont)/(totalCopy*(1-normCont)+totalCopyNormal*normCont)
	names(af)<-as.character(round(trueAF,3))
	return(af)
}
decontAF<-function(contAF,totalCopy,normCont=0,totalCopyNormal=2,type=c("mutation","SNPHet","SNPHomo")){
	type<-match.arg(type)
	totalCopy<-.checkLength(totalCopy,"totalCopy",contAF,"contAF")
	normCont<-.checkLength(normCont,"normCont",contAF,"contAF")
	totalCopyNormal<-.checkLength(totalCopyNormal,"totalCopyNormal",contAF,"contAF")
	if(any(normCont>1 | normCont<0)) stop("normCont must be between 0 and 1")
	const<-switch(type,"mutation"=0,"SNPHet"=1,SNPHomo=totalCopyNormal) #number of copies in the normal that contain the variant
	af<-(contAF*(totalCopy*(1-normCont)+2*normCont) -const*normCont)/(totalCopy*(1-normCont))
	names(af)<-as.character(round(contAF,3))
	return(af)
}

errorAF<-function(trueAF,seqError=0){
	seqError<-.checkLength(seqError,"seqError",trueAF,"trueAF")
	((3-4*seqError)*trueAF+seqError)/(3-2*seqError)
}

allAF<-function(totalCopy,normCont=0,totalCopyNormal=2,type=c("mutation","SNPHet","SNPHomo")){
    type<-match.arg(type)
		#     if(type=="mutation"){
		# possAlleles<-(1:totalCopy)/totalCopy
		# 	     pTumor<-sapply(normCont,function(p){totalCopy*(1-p)/(totalCopy*(1-p)+2*(p))})
		# AF<-lapply(pTumor,FUN=function(x){
		#             out<-possAlleles*x
		#             names(out)<-paste(1:totalCopy,"/",totalCopy,sep="")
		#             return(out)
		#         })
		#      }
    # if(type%in%c("SNPHet","SNPHomo") ){
        AF<-lapply(normCont,function(p){
            out<-sapply((0:totalCopy)/totalCopy, function(x){
				contAF(x,totalCopy,p,type=type,totalCopyNormal=totalCopyNormal)
			})
            names(out)<-paste(0:totalCopy,"/",totalCopy,sep="")
			if(type=="mutation") out<-out[-1]
            return(out)
        })
    # }
    names(AF)<-as.character(round(normCont,2))
    return(AF)

}

mleAF<-function(x,m,totalCopy,maxCopy=totalCopy,seqError=0,normCont=0){
	if(length(x)!=length(m)) stop("x and m must be the same length")
	if(length(normCont)>1 | length(seqError)>1 | length(totalCopy)>1 | length(maxCopy)>1){ 
		multiSample<-TRUE 	
		#note: this replicates everything, so a vector of same length. Messes up where not need to be a vector
		if(length(normCont)>1){if(length(x)!=length(normCont)) stop("x/m and normCont should be same length")} else normCont<-rep(normCont,length(x))
		if(length(seqError)>1){if(length(x)!=length(seqError)) stop("x/m and seqError should be same length")} else seqError<-rep(seqError,length(x))		
		if(length(totalCopy)>1){if(length(x)!=length(totalCopy)) stop("x/m and totalCopy should be same length")} else totalCopy<-rep(totalCopy,length(x))
		if(length(maxCopy)>1){if(length(x)!=length(maxCopy)) stop("x/m and maxCopy should be same length")} else maxCopy<-rep(maxCopy,length(x))		
	}
	else{multiSample<-FALSE}

	
	if(multiSample){
		r<-lapply(maxCopy,FUN=function(x,y){1:x} )
		alleles<-mapply(totalCopy,maxCopy,normCont,FUN=function(tot,mx,n){contAF(1:mx/tot,normCont=n,totalCopy=tot)},SIMPLIFY=FALSE)
		alleles<-mapply(seqError,alleles,FUN=function(e,a){errorAF(a,seqError=e)},SIMPLIFY=FALSE)
		out<-mapply(x,m,alleles,FUN=function(x,m,a){sapply(a,function(z){dbinom(x,size=m,prob=z)})},SIMPLIFY=FALSE)	
		out<-mapply(out,r,FUN=function(o,r){names(o)<-r; return(o)},SIMPLIFY=FALSE)
		assign<-mapply(r,sapply(out,which.max),FUN=function(r,wh){r[wh]})
	}
	else{
		r<-1:maxCopy 
		alleles<-contAF(r/totalCopy[1],normCont=normCont[1],totalCopy=totalCopy[1])
		alleles<-errorAF(alleles,seqError=seqError[1])
		out<-sapply(alleles,function(z){dbinom(x,size=m,prob=z)})
		colnames(out)<-r		
		assign<-r[apply(out,1,which.max)]
		assignAlleles<-factor(alleles[apply(out,1,which.max)],levels=alleles)
		afdf<-data.frame(tumorAF=r/totalCopy,AF=alleles,table(assignAlleles))
		row.names(afdf)<-NULL
		afdf<-afdf[,-3]
		
	}
	outList<-list(perLocationProb=out,assignments=data.frame(nCopies=assign,totalCopy=totalCopy,tumorAF=assign/totalCopy))
	if(!multiSample) outList<-c(outList,list(alleleSet=afdf))
	return(outList)
}