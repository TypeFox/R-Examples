plotMPPI <-
  function(x,threshold.model=0.01,threshold.variable=0.1,Figure=TRUE,cutoff=TRUE,useMC=FALSE){
    
    if(threshold.model<1) ListSelect <- which(x$BestModels$postProb>threshold.model) else{
      ListSelect <-1:threshold.model}
    
    nSelected <- length(ListSelect)
    
    
    if(!is.null(x$MAP.file)){
      SelectedModels<-list('Rank'=x$BestModels$Rank[ListSelect],'nVisits'=x$BestModels$nVisits[ListSelect],'ModeSize'=x$BestModels$ModeSize[ListSelect],'logCondPost'=x$BestModels$logCondPost[ListSelect],'jeffries'=x$BestModels$jeffries[ListSelect],'postProb'=x$BestModels$postProb[ListSelect],'modelName'=x$BestModels$modelName[ListSelect],'modelPosInX'=x$BestModels$modelPosInX[ListSelect])
    }else{
      SelectedModels<-list('Rank'=x$BestModels$Rank[ListSelect],'nVisits'=x$BestModels$nVisits[ListSelect],'ModeSize'=x$BestModels$ModeSize[ListSelect],'logCondPost'=x$BestModels$logCondPost[ListSelect],'jeffries'=x$BestModels$jeffries[ListSelect],'postProb'=x$BestModels$postProb[ListSelect],'modelName'=x$BestModels$modelName[ListSelect])
    }
    
    
    if(is.null(x$label.X)) label.X <- as.character(1:x$p)
    if(is.null(x$label.Y)) Pheno <- paste("Y",1:x$q,sep="",collapse="_") else Pheno <- paste("Y: ",paste(x$label.Y,collapse="/"),sep="") 
    
    if (useMC) {
        NameMarg <- file.path(x$path.output, paste(x$root.file.output,"output_marg_prob_incl_mc.txt",sep="_"))
    } else {
        NameMarg <- file.path(x$path.output, paste(x$root.file.output,"output_marg_prob_incl.txt",sep="_"))
    }
    Marg <- read.table(NameMarg,header=TRUE)
    
    lab.X <- "predictor"
    if(Figure==TRUE){
      
     # x11(width=13,height=6)
       
      plot(Marg$Predictor,type='h',Marg$Marg_Prob_Incl,xlab=lab.X,ylab='MPPI: Marginal Posterior Probability Inclusion',col='blue',xaxt='n',lty=1,main=Pheno,cex.main=1.5,axes=FALSE)
      atx <- seq(0,1,by=0.2)
      axis(2,at=atx,labels=atx)
      
      if(!is.null(x$MAP.file)){
        if(is.data.frame(x$MAP.file)){SNPLabels <- x$MAP.file}else{
          NameMap.file <- file.path(x$path.input,x$MAP.file)
          SNPLabels <- read.table(NameMap.file,header=TRUE,stringsAsFactors = FALSE)}
        label.X <- SNPLabels$SNPName
        FullSNPList <- as.numeric(unique(unlist(strsplit(SelectedModels$modelPosInX," "))))
        UniqueSNPS<-list('PosInX'=FullSNPList,'chr'=SNPLabels$Chr[FullSNPList],'Posn'=SNPLabels$Posn[FullSNPList],'rsName'=SNPLabels$SNPName[FullSNPList],'margProb'=Marg$Marg_Prob_Incl[FullSNPList])
        Order <- sort(UniqueSNPS$margProb,index.return=T,decreasing=T)$ix
        
        SortedUniqueSNPS<-list('PosInX'=UniqueSNPS$PosInX[Order],'chr'=UniqueSNPS$chr[Order],'Posn'=UniqueSNPS$Posn[Order],'rsName'=UniqueSNPS$rsName[Order],'margProb'=UniqueSNPS$margProb[Order])
        
        ListChr <- unique(SNPLabels$Chr)
        nChr <- length(ListChr)
        
        ColFill <- rep(c('lightgray','lightblue'),(floor(nChr/2)+1))
        AtX <- c()
        for(CHR in 1:nChr){
          LowerBd <- min(which(SNPLabels$Chr==CHR))
          UpperBd <- max(which(SNPLabels$Chr==CHR))            
          labTxt=paste('Chr ',ListChr[CHR])
          text(x=LowerBd+(UpperBd-LowerBd)*0.5,y=1.1,labels=labTxt,pos=4,cex=2)
          rect(LowerBd,0,UpperBd,1,col=ColFill[CHR],border=F)
          AtX[CHR] <- LowerBd+(UpperBd-LowerBd)*0.5
        }
        
        abline(h=0,col='blue')
        TheshPoints <- which(Marg$Marg_Prob_Incl>0.0001)
        lines(Marg$Predictor[TheshPoints],Marg$Marg_Prob_Incl[TheshPoints],col='blue',type='h')
        Highlight <- which(SortedUniqueSNPS$margProb>threshold.variable)
        LabelSnps <- SortedUniqueSNPS$rsName[Highlight]
        PosSNPs <- SortedUniqueSNPS$PosInX[Highlight]
        PbtySNPs <- SortedUniqueSNPS$margProb[Highlight]
        points(PosSNPs,PbtySNPs,col='red',pch=19)
        text(x=PosSNPs,y=PbtySNPs,labels=LabelSnps,srt=90,pos=2,cex=1.3)
        XLAb <- paste("Chr",1:nChr)
        axis(1,at=AtX,labels=XLAb,srt=90,cex=1.3)
        Highlight.MPI <- which(Marg$Marg_Prob_Incl>threshold.variable)
     
      }else{
          if(x$p>20){
            position <- seq(1,x$p,by=round(x$p/20))
          XLAb <- paste(" ",position)
          axis(1,at=position,labels=XLAb,srt=180,cex=1.3,las=2) 
         }
        }
      
      x1=par('usr')[1]
      x2=par('usr')[2]
      y1=par('usr')[3]
      y2=par('usr')[4]
      xleg=x1+(x2-x1)*0.0
      yleg=y2+(y2-y1)*0.8
      
      if(cutoff==TRUE){ abline(h=threshold.variable,col='black',lty=2)
      
      
      legend("topleft", c("In best models", "MPPI cut-off"), col = c("red","black"),
             text.col = c("red","black"), lty = c(NA,2), pch = c(19, NA),
             merge = TRUE, bg = 'gray90')#,inset=c(0,-0.2))
      
      }else{
      
      
      legend("topleft", c("In best models"), col = c("red"),
             text.col = c("red"), pch = c(19),
              bg = 'gray90')#,inset=c(0,-0.2))
      }
      grid(nx=0,ny=6,lty=3,col="gray")
      
    #
    }
    
    Highlight <- which(Marg$Marg_Prob_Incl>threshold.variable)
    
    if(length(Highlight)==0){
      cat("none variable has a MPPI greater than ",threshold.variable,"\n")
      return(c(SelectedModels,list(var.TOP.MPI=NULL,var.MPI=NULL)))
    }else{ 
      if((Figure==TRUE) & is.null(x$MAP.file)){
        
        
        PosSNPs <- Highlight
        PbtySNPs <- Marg$Marg_Prob_Incl[Highlight]
        text(x=Highlight,y=PbtySNPs,labels=label.X[Highlight],pos=1,cex=1.0,srt=90,offset=0.9,adj=c(1,1))
      }
      ##variable MPI > threshold.variable and in best model > threshold.model --> highlight in red
      varX <- as.character(1:x$p)
      
      if(!is.null(x$MAP.file)){
        if(is.data.frame(x$MAP.file)){SNPLabels <- x$MAP.file}else{
          NameMap.file <- file.path(x$path.input,x$MAP.file)
          SNPLabels <- read.table(NameMap.file,header=TRUE,stringsAsFactors = FALSE)}
        varX <- label.X <- SNPLabels$SNPName}
      
      
      varHighlight <- which(Marg$Marg_Prob_Incl>threshold.variable&(varX%in%unique(unlist(strsplit(SelectedModels$modelName," ")))))
      
      
      if(length(varHighlight)>0&(Figure==TRUE)&is.null(x$MAP.file)){
        PosSNPs <- varHighlight
        PbtySNPs <- Marg$Marg_Prob_Incl[varHighlight]
        points(PosSNPs,PbtySNPs,col='red',pch=19)
      }
      
      return(c(SelectedModels,list(var.TOP.MPI=as.character(label.X[varHighlight]),var.MPI=as.character(label.X[Highlight]))))
    }
    
  }

#dev.off()
