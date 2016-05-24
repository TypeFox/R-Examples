PlotCGH <- function(array,chrom=NULL,nuc=NULL,main="", scaleEachChrom=TRUE){
  if(is.null(chrom)){
    chrom <- rep(1,length(array))
    warning("Since chrom was not entered, PlotCGH assumed that all CGH spots in array are on the same chromosome.")
  }
#  if(!is.numeric(chrom)) stop("Chrom must be numeric.")
  if(is.null(nuc)){
    nuc <- rep(NA, length(chrom))
    for(i in unique(chrom)){
      nuc[chrom==i] <- 1:sum(chrom==i)
    }
  }
  scaledarray <- numeric(length(array))
  if(scaleEachChrom){
    for(i in (unique(chrom))) scaledarray[chrom==i] <- array[chrom==i]/(1.1*max(abs(array[chrom==i])))
  } else {
    scaledarray <- array/(.9*max(abs(array)))
  }
  plot.CGH.FL.Single(scaledarray,chrom,nuc,main)
}

plot.CGH.FL.Single<-function(array, chr, nucposi, main=""){
  if(length(array)!=length(chr) || length(array)!=length(nucposi)) stop("Array, chrom, and nuc must all have the same length (or chrom & nuc can be NULL).")
  plot(0,0,type="n",axes=F,ylim=c(0,length(unique(chr))+1),xlim=c(-.05*max(nucposi), max(nucposi)),ylab="",xlab="",main=main,cex.main=1)
  for(j in 1:length((unique(chr)))){
    chrj <- (unique(chr))[j]
    jp=length((unique(chr)))-j+1
    nuc=nucposi[chr==chrj]
    y=array[chr==chrj]
    y[is.na(y)]<-0
    yposi=ynega=y
    yposi[y<0]<-0
    ynega[y>0]<-0
    pick<-(1:length(y))[y!=0]
    if(length(pick)>0){
      segments(nuc[pick],jp,nuc[pick],jp+yposi[pick],col=2)
      segments(nuc[pick],jp,nuc[pick],jp+ynega[pick],col=3)
    }
    segments(0,jp,max(nuc),jp)
    text(-.05*max(nucposi),jp,labels=chrj,cex=.7)
  }
}
