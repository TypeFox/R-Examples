## Julia Bischof
## 10-09-2015

#library(gplots)

plotGeneComb<-function(geneComb.tab=NULL,color=c("gray97","darkblue"), withNA=TRUE,title=NULL, PDF=NULL,...){
  if(length(geneComb.tab)==0){
    stop('--> Gene combination table is missing')
  }
  if(length(color)==1){
    color<-c("white",color)
  }
  if(max(geneComb.tab,na.rm=T)>1){
    geneComb.tab<-geneComb.tab/sum(geneComb.tab,na.rm=T)
  }
  
  if(nrow(geneComb.tab)<ncol(geneComb.tab)){
    geneComb.tab<-(t(geneComb.tab))
  }
  if(length(PDF)>0){
    pdf(file = paste(PDF,"_",substr(colnames(geneComb.tab)[1],1,4),"-",substr(rownames(geneComb.tab)[1],1,4),"-combinations.pdf",sep=""),
        width = nrow(geneComb.tab)+8,
        height = ncol(geneComb.tab)+3,
        pointsize = if(sqrt(min(dim(geneComb.tab))*min(dim(geneComb.tab)))<18){18}else{round(sqrt(min(dim(geneComb.tab))*min(dim(geneComb.tab)))+10,-1)})
  }
  if(withNA==FALSE){
    geneComb.tab<-geneComb.tab[-which(rownames(geneComb.tab)=="NA"),-which(colnames(geneComb.tab)=="NA")]
  }
  par(oma=c(nchar(colnames(geneComb.tab)[1])*0.6,1,2,nchar(colnames(geneComb.tab)[1])*0.6))
  heatmap.2(t(as.matrix(geneComb.tab)),col=colorRampPalette(c(color[1],color[2])),Rowv=T,Colv=T,key=T,density.info="none",
            main=if(length(title)==0){paste(substr(colnames(geneComb.tab)[1],1,4), " & ", substr(rownames(geneComb.tab)[1],1,4)," combinations",sep="")}else{title}, 
            trace="none",cexCol=1.3,cexRow=1.3,dendrogram="none",key.xlab="proportion",key.par=list(cex=0.9),cex.main=1.4)
  
  if(length(PDF)>0){
    dev.off()
  }
}


