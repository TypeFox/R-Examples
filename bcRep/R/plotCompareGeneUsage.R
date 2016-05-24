## Julia Bischof
## 15-10-2015

plotCompareGeneUsage<-function(comp.tab=NULL, color=c("gray97","darkblue"), title=NULL, PDF=NULL){
  if(length(comp.tab)==0){
    stop("--> Comparison table is missing")
  }
  if(length(color)==1){
    color<-c("white",color)
  }
  if(length(PDF)>0){
    pdf(paste(PDF,"_Comparison_Gene-usage.pdf",sep=""),
        width=if(ncol(comp.tab)>10){ncol(comp.tab)*0.6}else{12},
        height=if(ncol(comp.tab)>10){ncol(comp.tab)/3}else{10},
        pointsize=if(ncol(comp.tab)>10 || nrow(comp.tab)>10){28}else{16})
  }
  par(oma=c(4*nrow(comp.tab)/5,0,0,4*nrow(comp.tab)/5))
  heatmap.2(as.matrix(comp.tab),col=colorRampPalette(c(color[1],color[2])),Rowv=T,Colv=T,key=T,density.info="none",
            main=if(length(title)==0){paste(substr(colnames(comp.tab)[1],1,4)," usage",sep="")}else{title}, 
            trace="none",dendrogram="both",key.xlab=if(max(comp.tab,na.rm=T)>1){"quantity"}else{"proportion"},
            key.par=list(cex=0.9), cexCol=1.1, cexRow=1.1)
  if(length(PDF)>0){
    dev.off()
  }
}


