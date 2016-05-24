#'Manhattan Plot for meta-analysis results 
#'
#'\code{mhplot} returns a pdf file containing the Manhattan plot of meta-analysis results.
#'
#'The function is useful for a quick overview of meta-analysis results.
#'Two different types of input are allowed: 
#'\itemize{ 
#'\item a data frame containing chromosome, position and p-values to be plotted.  
#'\item a file name for retrieving information from file. Note that separator is set to white space,
#'can be changed accordingly if needed. Special character "*" is allowed for selecting more than one file (e.g. in case there is
#'one file for each chromosome).
#'}
#'
#'
#'@param res.file File name for the file(s) containing meta-analysis results. For multiple files (e.g. one file for each chromosome)
#' the special character "*" is allowed. This is an alternative to \code{data}.
#'@param data Dataset containing meta-analysis results. This is an alternative to \code{file}, if results are already loaded
#' in the R workspace.
#'@param col Choice of two colours for the plot. They will be alternating for different chromosomes.
#'@param out.file File name for the output.
#'@param main Title to appear on plot.   
#'@param plab p-values label i.e. column name for the p-values column.
#'@param CHRlab Chromosome label i.e. column name for the chromosome numbers column.
#'@param POSlab Position label i.e. column name for the position numbers column.
#'@param sig.p Significant p-values level. 
#'@param sug.p Suggestive p-values level.
#'@param sep Separator used in \code{file}. Default is white space. 
#'
#'@return The output is a pdf file containing the Manhattan Plot.
#'@import mvtnorm expm 
#'@export


mhplot <- function(res.file=NULL,data=NULL,col=c("darkgrey","lightgrey"),out.file="mhplot.pdf",main="MHplot",plab="p_value",CHRlab="chr",POSlab="Position",
             sig.p=5e-8,sug.p=5e-7,sep=" "){

  if(length(res.file)!=0){
  files=system(paste("ls ", res.file, sep=""), intern=T)
  res=c()    #input data
  if(length(files)>1){
    for(i in 1:length(files)){
      res.he=scan(files[i], nlines=1, "char", quiet=T)
      he.ind1=grep(CHRlab, res.he)
      he.ind2=grep(POSlab, res.he)
      he.ind3=grep(plab, res.he)
      res.tmp=read.table(pipe(paste("cut -d \"",sep,"\" -f ",he.ind1,",",he.ind2,",",he.ind3," ",files[i],sep="")), header=T, stringsAsFactors=F)
    plim=10000/dim(res.tmp)[1]
    if(length(which(res.tmp[,plab]>plim))>10000){
      idx=c(which(res.tmp[,plab]<plim),sample(which(res.tmp[,plab]>plim),size=10000))  #lighten the burden of non-significant points
      res.tmp=res.tmp[idx,]}
      res=rbind(res, res.tmp)
  }}  
  if(length(files)==1){
    res=read.table(res.file,header=T,stringsAsFactors=F,sep=sep)
    plim=200000/dim(res)[1]
    if(length(which(res[,plab]>plim))>200000){
      idx=c(which(res[,plab]<plim),sample(which(res[,plab]>plim),size=200000))
      res=res[idx,]
    }
  }
  }else{
    res=data[,c(CHRlab, POSlab, plab)]
  }

    chromosomes=unique(res[,CHRlab])
    if(any(chromosomes=="X")){
       res[which(res[,CHRlab]=="X"),CHRlab]=23
      chromosomes[chromosomes=="X"]=23
      chromosomes=as.numeric(chromosomes)
          }
    chromosomes=sort(chromosomes)
    res=res[which(!is.na(res[,POSlab])),]
    maxs=by(res[,POSlab],res[,CHRlab],max)
    nms=names(maxs)
    maxs=as.vector(maxs)
    names(maxs)=nms
    offset=cumsum(as.numeric(maxs))
    if(max(-log10(res[,plab]),na.rm=T)>10){
      ylim=c(0,max(-log10(res[,plab]),na.rm=T))
    }else{
      ylim=c(0,10)}
    offset=c(0,offset)
    pdf(paste(out.file,".pdf",sep=""),width=21)
    plot(1,xlim=c(0,max(offset)),ylim=ylim,type="n",xaxt="n",ylab="-log10(p)", xlab="Chromosome", main=main)
    for(i in 1:length(chromosomes)){
      res.tmp=res[which(res[,CHRlab]==chromosomes[i]),]
      points(x=res.tmp[,POSlab]+offset[i],y=-log10(res.tmp[,plab]),pch=20,col=col[(i%%2)+1])
      res.sig=res.tmp[which(res.tmp[,plab]<=sig.p),]
      if(nrow(res.sig)>0){
        points(x=res.sig[,POSlab]+offset[i],y=-log10(res.sig[,plab]),pch=20,col="red")
      }
      res.sug=res.tmp[which(res.tmp[,plab]<=sug.p & res.tmp[,plab]>sig.p),]
      if(nrow(res.sug)>0){
        points(x=res.sug[,POSlab]+offset[i],y=-log10(res.sug[,plab]),pch=20,col="blue")
      }
    }
    abline(h=-log10(sig.p),col="red")
    abline(h=-log10(sug.p),col="blue")
    axis(side=1,at=rowMeans(cbind(offset[-length(offset)],offset[-1])),labels=chromosomes)
    dev.off()
  }
