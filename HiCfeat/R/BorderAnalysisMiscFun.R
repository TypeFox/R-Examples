# Raphael Mourad
# Cuvier team, LBME lab
# 15/04/2015


# Read bed file containing ChIP-seq peak data
readGFBed<- function(GFBedFile,seqInfoChr){

 if(class(seqInfoChr)!="Seqinfo"){print("seqInfoChr is not a seqinfo object!");return(0)}
 if(length(grep(pattern='.bed',x=GFBedFile))==0){print("GFBedFile is not a bed file!");return(0)}

 dataGF=read.table(GFBedFile,sep='\t',header=F)[,c(1:3)]
 dataGF2=dataGF[dataGF[,1]%in%seqnames(seqInfoChr),]
 chr=as.character(dataGF2[,1])
 posLeft=as.numeric(dataGF2[,2])
 posRight=as.numeric(dataGF2[,3])

 GenomicFeature.GR=NULL
 for(i in 1:length(seqnames(seqInfoChr))){
  chri=seqnames(seqInfoChr)[i]
  if(sum(chr==chri)){
   GenomicFeature.GRi <- GRanges(seqnames=chri,IRanges(start=posLeft[chr==chri],end=posRight[chr==chri]))
   if(i==1){
    GenomicFeature.GR=GenomicFeature.GRi
   }else{
    GenomicFeature.GR=c(GenomicFeature.GR,GenomicFeature.GRi)
   }
  }
 }

 seqlevels(GenomicFeature.GR) = seqlevels(seqInfoChr)
 seqinfo(GenomicFeature.GR)=seqInfoChr

 #print("Bed file read")

 return(GenomicFeature.GR)
}

# Read wig file containing ChIP-seq peak data
readGFWig<- function(GFWigFile,seqInfoChr){

 if(class(seqInfoChr)!="Seqinfo"){print("seqInfoChr is not a seqinfo object!");return(0)}
 if(length(grep(pattern='.wig',x=GFWigFile))==0 & length(grep(pattern='.bw',x=GFWigFile))==0){print("GFWigFile is not a wig or bigwig file!");return(0)}

 dataGF.GR=import(GFWigFile,asRangedData=FALSE)
 dataGF=as.data.frame(dataGF.GR)
 dataGF2=dataGF[dataGF[,1]%in%seqnames(seqInfoChr),]
 chr=as.character(dataGF2[,1])
 posLeft=as.numeric(dataGF2[,2])
 posRight=as.numeric(dataGF2[,3])
 score=as.numeric(dataGF2[,6])
 #GenomicFeature.GR <- GRanges(seqnames=chr,IRanges(start=posLeft,end=posRight),score=score)

 GenomicFeature.GR=NULL
 for(i in 1:length(seqnames(seqInfoChr))){
  chri=seqnames(seqInfoChr)[i]
  if(sum(chr==chri)){
   GenomicFeature.GRi <- GRanges(seqnames=chri,IRanges(start=posLeft[chr==chri],end=posRight[chr==chri]),score=score[chr==chri])
   if(i==1){
    GenomicFeature.GR=GenomicFeature.GRi
   }else{
    GenomicFeature.GR=c(GenomicFeature.GR,GenomicFeature.GRi)
   }
  }
 }

 seqlevels(GenomicFeature.GR) = seqlevels(seqInfoChr)
 seqinfo(GenomicFeature.GR)=seqInfoChr

 #print("Bed file read")

 return(GenomicFeature.GR)
}



# Read bed file containing domain data
readDomBed<- function(domainBedFile,seqInfoChr){

 if(class(seqInfoChr)!="Seqinfo"){print("seqInfoChr is not a seqinfo object!");return(0)}
 if(length(grep(pattern='.bed',x=domainBedFile))==0){print("GFBedFile is not a bed file!");return(0)}

 dataDomains <-read.table(domainBedFile,sep='\t',header=F)
 dataDomains2=dataDomains[dataDomains[,1]%in%seqnames(seqInfoChr),]
 chr=as.character(dataDomains2[,1])
 posLeft=as.numeric(dataDomains2[,2])
 posRight=as.numeric(dataDomains2[,3])

 Domains.GR=GRanges(seqnames=chr,IRanges(start=posLeft,end=posRight),seqinfo=seqInfoChr)


 Domains.GR=NULL
 for(i in 1:length(seqnames(seqInfoChr))){
  chri=seqnames(seqInfoChr)[i]
  if(sum(chr==chri)){
   Domains.GRi <- GRanges(seqnames=chri,IRanges(start=posLeft[chr==chri],end=posRight[chr==chri]))
   if(i==1){
    Domains.GR=Domains.GRi
   }else{
    Domains.GR=c(Domains.GR,Domains.GRi)
   }
  }
 }


 seqlevels(Domains.GR) = seqlevels(seqInfoChr)
 seqinfo(Domains.GR)=seqInfoChr

 if(ncol(dataDomains)==4){
  Domains.GR$DomainType=dataDomains[,4]
 }

 #print("Bed file read")

 return(Domains.GR)
}


# Average values by binning of a GRanges object
averagePerBin <- function(x, binsize, mcolname)
{
 if(!is(x, "GenomicRanges")){stop("'x' must be a GenomicRanges object")}
 if(any(is.na(seqlengths(x)))){stop("'seqlengths(x)' contains NAs")}

 bins <- IRangesList(lapply(seqlengths(x),function(seqlen)IRanges(breakInChunks(seqlen, binsize))))

 cvg <- coverage(x, weight=mcolname)

 views_list <- RleViewsList(lapply(names(cvg),function(seqname)Views(cvg[[seqname]], bins[[seqname]])))

 averageBin=NULL
 for(i in 1:length(views_list)){
  averageBin=c(averageBin,viewMeans(views_list)[[i]])
 }

 return(averageBin)
}


