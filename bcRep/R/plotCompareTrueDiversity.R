## Julia Bischof
## 10-15-2015

plotCompareTrueDiversity<-function(comp.tab=NULL, mean.plot=T, colors=NULL, title=NULL, PDF=NULL){
  if(length(comp.tab)==0){
    stop("--> Comparison table is missing")
  }
  if(length(title)>0){
    plottitle<-title
  }else{
    plottitle<-NULL
  }
  
  uniseqlength<-vector()
  for(i in 2:length(comp.tab)){
    uniseqlength<-c(uniseqlength,names(comp.tab[[i]]))
  }
  uniseqlength<-unique(uniseqlength)
  orderCDR3<-as.numeric(unlist(apply(data.frame(uniseqlength),1,function(x){strsplit(x, split=" = ")[[1]][2]})))
  uniseqlength<-uniseqlength[order(orderCDR3)]
  
  if(length(colors)!=(length(comp.tab)-1)){
    color<-rainbow(n=(length(comp.tab)-1),start=0.2, end=0.8)
  }else{
    color<-colors
  }
  
  nrleg<-ceiling((length(comp.tab)-1)/4)
  
  
  if(mean.plot==T){
    if(length(PDF)>0){
      pdf(paste(PDF,"_Comparison_True-diversity_q",as.numeric(comp.tab$True_diversity_order),".pdf",sep=""),
        width=9,height=7, pointsize =16)
    }
    par(mfrow=c(1,1),mar=c(5,5,2,2), oma=c(1,1,4,1))    

    maxy<-round(max(unlist(comp.tab),na.rm=T),0)+1
    plot(1,1, col="white",xlab="CDR3 length [AA]", ylab="Mean diversity", xlim=c(1,length(uniseqlength)+5), 
         ylim=c(0, maxy))
    for(i in 2:length(comp.tab)){
      m<-vector()
      s<-vector()
      for(j in 1:length(comp.tab[[i]])){
        m<-c(m, mean(as.numeric(comp.tab[[i]][[j]]), na.rm=T))
        s<-c(s, sd(as.numeric(comp.tab[[i]][[j]]), na.rm=T))
      }
      points(x = as.numeric(apply(data.frame(names(comp.tab[[i]])),1,function(x){strsplit(x, split="= ")[[1]][2]})), m, pch=16, col=color[i-1], cex=1.3)
      lines(x = as.numeric(apply(data.frame(names(comp.tab[[i]])),1,function(x){strsplit(x, split="= ")[[1]][2]})), m, lwd=2, col=color[i-1])
      CI.up = as.numeric(m)+as.numeric(s)
      CI.dn = as.numeric(m)-as.numeric(s)
      options(warn=-1)
      arrows(as.numeric(apply(data.frame(names(comp.tab[[i]])),1,function(x){strsplit(x, split="= ")[[1]][2]})),CI.dn,
             as.numeric(apply(data.frame(names(comp.tab[[i]])),1,function(x){strsplit(x, split="= ")[[1]][2]})),CI.up,code=3,length=0.05, angle=90,col=color[i-1])
      options(warn=0)
    }
    newnames<-names(comp.tab)[2:length(comp.tab)]
    if(nrleg>1){
      for(j in 1:nrleg){
        legend("topright",col=color[((j*4)-3):(j*4)],newnames[((j*4)-3):(j*4)],lty=1, lwd=3, cex=1.1, y.intersp = 1.2)
      }
    }else{
      legend("topright",col=color,newnames,lty=1, lwd=3, cex=1, y.intersp = 1.1)
    }    
  }else{
    if(length(PDF)>0){
      pdf(file = paste(PDF,"_Comparison_True-diversity_q",as.numeric(comp.tab$True_diversity_order),".pdf",sep=""),
          width = sqrt((length(uniseqlength)+nrleg+1))*2.5,
          height = sqrt((length(uniseqlength)+nrleg+1))*2.3, 
          pointsize = if(sqrt(length(uniseqlength)+nrleg)*2.1<12){12}else{sqrt(length(uniseqlength)+nrleg)*2.1})
    }
    par(mfrow=c(ceiling(sqrt((length(uniseqlength)+nrleg+1))),ceiling(sqrt((length(uniseqlength)+nrleg+1)))),
        mar=c(3,4,7,2),oma=c(2,2,10,3))
    maxy<-round(max(unlist(comp.tab),na.rm=T),0)+1
    for(i in 1:length(uniseqlength)){
      plot(1,1,col="white",xlab=" ",xaxt="n",main=uniseqlength[i],ylim=c(0,maxy),cex=0.8, ylab=" ", 
           xlim=c(0.8,as.numeric(strsplit(uniseqlength[i],split=" = ")[[1]][2])+0.2),cex.axis=0.9)
      if(as.numeric(strsplit(uniseqlength[i],split=" = ")[[1]][2])<=5){
        axis(1,at=seq(1,as.numeric(strsplit(uniseqlength[i],split=" = ")[[1]][2]),1),seq(1,as.numeric(strsplit(uniseqlength[i],split=" = ")[[1]][2]),1),cex.axis=0.9)
      }else if(as.numeric(strsplit(uniseqlength[i],split=" = ")[[1]][2])<=14){
        axis(1,at=seq(1,as.numeric(strsplit(uniseqlength[i],split=" = ")[[1]][2]),2),seq(1,as.numeric(strsplit(uniseqlength[i],split=" = ")[[1]][2]),2),cex.axis=0.9)
      }else{
        axis(1,at=c(1,round(as.numeric(strsplit(uniseqlength[i],split=" = ")[[1]][2])*seq(0.25,1,0.25),0)),c(1,round(as.numeric(strsplit(uniseqlength[i],split=" = ")[[1]][2])*seq(0.25,1,0.25),0)),cex.axis=0.9)
      }
      mtext(side = 2, at = maxy/2, las=3, line = 2, cex=0.7,
            text = if(comp.tab$True_diversity_order==0){"Richness"}else{"Diversity"})
      mtext(side = 1, line = 2, cex=0.7, text = "Position")
      for(j in 2:length(comp.tab)){
        if(as.numeric(strsplit(uniseqlength[i],split=" = ")[[1]][2])==1){
          if(uniseqlength[i] %in% names(comp.tab[[j]])){
            points(as.numeric(comp.tab[[j]][[which(names(comp.tab[[j]])==uniseqlength[i])]]),
                   col=color[j-1], pch=20)
          }
        }
        if(uniseqlength[i] %in% names(comp.tab[[j]])){
          lines(as.numeric(comp.tab[[j]][[which(names(comp.tab[[j]])==uniseqlength[i])]]),
                col=color[j-1], lwd=2)
        }
      }
    }
    newnames<-names(comp.tab)[2:length(comp.tab)]
    if(nrleg>1){
      for(j in 1:nrleg){
        plot(1,1,col="white",axes=F,xlab="",ylab="")
        legend("center",col=color[((j*4)-3):(j*4)],newnames[((j*4)-3):(j*4)],lty=1, lwd=3, cex=1.1, y.intersp = 1.2, box.col = "white")
      }
    }else{
      plot(1,1,col="white",axes=F,xlab="",ylab="")
      legend("center",col=color,newnames,lty=1, lwd=3, cex=1, y.intersp = 1.1, box.col = "white")
    }
  }
  title(if(length(plottitle)>0){plottitle}else{paste("True diversity, q = ", as.numeric(comp.tab$True_diversity_order),sep="")},outer=T, cex.main=2)
  if(length(PDF)>0){
    dev.off()
  }
}
