plotSPs<-function(dm, sampleID=NA,cex=0.5, legend="CN_Estimate", orderBy="chr", rawAF=F){
  
  dm[,"%maxP"]=dm[,"%maxP"]-min(dm[,"%maxP"],na.rm=T)+1;
  dm[,"%maxP"]=50*dm[,"%maxP"]/max(dm[,"%maxP"],na.rm=T);
  
  keep=which(!is.na(dm[,"SP"]),); dm=dm[keep,];
  ia=order(dm[,"startpos"]);dm=dm[ia,];
  ia=order(dm[,orderBy]);dm=dm[ia,];
  ia=order(dm[,"SP"],decreasing = TRUE);dm=dm[ia,];
  
  maxPloidy=max(dm[,"PM"],na.rm=TRUE);
  yticklab=c(1:maxPloidy,seq(0,1,by=0.1));
  at=c(sort(c(1:maxPloidy)*-0.1),seq(0,1,by=0.1));
  
  dm=.addColumn(dm,"AF_Tumor_Adjusted",NA);
  if(!rawAF){
    iEq2=which(dm[,"PM_B"]==dm[,"PN_B"])
    iEq3=setdiff(1:nrow(dm),iEq2);
    if(!isempty(iEq3)){ ##Equation 3 informative
      dm[iEq3,"AF_Tumor_Adjusted"]=(dm[iEq3,"AF_Tumor"]*dm[iEq3,"CN_Estimate"]-dm[iEq3,"PN_B"])/(dm[iEq3,"PM_B"]-dm[iEq3,"PN_B"])
      dm[iEq3,"AF_Tumor_Adjusted"]=dm[iEq3,"AF_Tumor_Adjusted"]*(dm[iEq3,"PM_B"]/dm[iEq3,"PM"])
    }
    if(!isempty(iEq2)){##Equation 3 not informative --> use equation 2
      dm[iEq2,"AF_Tumor_Adjusted"]=(dm[iEq2,"CN_Estimate"]-2)/(dm[iEq2,"PM"]-2)
    }
    adjusted="Adjusted"
  }else{
    dm[,"AF_Tumor_Adjusted"]=dm[,"AF_Tumor"]
    adjusted=""
  }
  
  
  par(xpd=T, cex=cex, cex.axis=1/cex,cex.lab=1/cex, cex.main=1/cex,mar=par()$mar+c(0,0.5,0,4.2));
  #   plot.new(); 
  par(xpd=FALSE)
  
  plot(c(1:length(at)),at, col="white",pch=8,xlim=c(0,nrow(dm)),
       yaxt="n", bty="L", main=sampleID, xlab="Mutation", 
       ylab=paste("Copy-number <--->",adjusted,"Allele-frequency and SP size"));
  axis(2, at, yticklab) 
  
  legend1=.plotSPPerVar(dm,17,1,0,legend);
  
  x=fliplr(gray.colors(50))
  # norm=1/length(x);
  for (k in 1:nrow(dm)){
    ci=max(1,ceil(dm[k,"%maxP"])); #/norm
    matpoints(k,dm[k,"SP"],pch=15,col=x[ci]);
    if (k==1){
      legend1$text[length(legend1$text)+1]="SP";
      legend1$col[length(legend1$text)]=x[ci];
      par(xpd=TRUE)
      legend("topright",legend1$text,fill=legend1$col,inset=c(-0.1,-0.02),cex=0.75/cex,bty = "n")
    }
  }
  .plotSPPerVar(dm,8,0,0,legend);
  .plotSPPerVar(dm,8,0,1,legend);
  .plotSPPerVar(dm,8,1,1,legend);
  lines(c(0,nrow(dm)),c(0,0),col="black");
  
  #ylim([-0.12*maxPloidy,1.1]);
}


.plotSPPerVar<-function(dm,lineType,lohFlag,cnFlag,var){
  if(var=="CN_Estimate"){
    .plotSPPerCopyNumber(dm,lineType,lohFlag,cnFlag)
  }else if(var=="chr"){
    .plotSPPerChr(dm,lineType,lohFlag,cnFlag)
  }
}

.plotSPPerChr<-function(dm,lineType,lohFlag,cnFlag){
  x=rainbow(40);
  legend1=list("text"=c(),"col"=c());
  maxploidy=max(dm[,"PM"],na.rm=TRUE)+1;
  for (i in 1:22){
    idx=which(dm[,"chr"]==i & ((lohFlag & dm[,"PN_B"]==1) |(!lohFlag & dm[,"PN_B"]==0)) );
    if (!isempty(idx)){
      if (!cnFlag){
        matpoints(idx,dm[idx,"AF_Tumor_Adjusted"],pch=lineType,col=x[i]);
      }else{
        if (any(!is.na(dm[idx,"PM"]))){
          matpoints(idx,0.1*(dm[idx,"PM"]-maxploidy),pch=20,col=x[i]);
          idy=intersect(idx,which(dm[,"PM"]!=dm[,"PM_cnv"]))
          if(!isempty(idy)){
            matpoints(idy,0.1*(dm[idy,"PM_cnv"]-maxploidy),pch=3,col=x[i]);
          }
        }
      }
    }
    legend1$text[i]=paste("chr",i);
    legend1$col[i]=x[i];
  }
  return(legend1);
}


.plotSPPerCopyNumber<-function(dm,lineType,lohFlag,cnFlag){
  x=c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C", "#FDBF6F", "#FF7F00");# "#CAB2D6");
  legend1=list("text"=c(),"col"=c());
  maxploidy=max(dm[,"PM"],na.rm=TRUE)+1;
  for (i in 1:length(x)){
    idx=which(round(dm[,"CN_Estimate"])==i & ((lohFlag & dm[,"PN_B"]==1) |(!lohFlag & dm[,"PN_B"]==0)) );
    if (!isempty(idx)){
      if (!cnFlag){
        matpoints(idx,dm[idx,"AF_Tumor_Adjusted"],pch=lineType,col=x[i]);
      }else{
        if (any(!is.na(dm[idx,"PM"]))){
          matpoints(idx,0.1*(dm[idx,"PM"]-maxploidy),pch=20,col=x[i]);
          idy=intersect(idx,which(dm[,"PM"]!=dm[,"PM_cnv"]))
          if(!isempty(idy)){
            matpoints(idy,0.1*(dm[idy,"PM_cnv"]-maxploidy),pch=3,col=x[i]);
          }
        }
      }
    }
    legend1$text[i]=paste("CN",i);
    legend1$col[i]=x[i];
  }
  return(legend1);
}
