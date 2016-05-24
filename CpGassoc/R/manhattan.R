manhattan <-
function(x,cpgname,chr,pos,save.plot=NULL,file.type="pdf",popup.pdf=FALSE,eps.size=c(15,5),main.title=NULL,cpg.labels=NULL,chr.list=NULL,color.list=NULL,point.size=NULL,...) {  
    if(!is.null(save.plot)){
      if(!(file.type %in% c("eps","pdf"))) {
             stop("Incorrect file type. Must be pdf or eps\n")
               }
      if(file.type=="eps"){
        postscript(paste(save.plot[1],".eps",sep=""),horizontal = FALSE, 
                onefile = FALSE, paper = "special",width=eps.size[1],height=eps.size[2])
             }
      if(file.type=="pdf" & !popup.pdf) {
         pdf(file=paste(save.plot[1],".pdf",sep=""))
          }
        }
    test<-x$results[,c(1,3)]
    chr<-gsub("[[:space:]]","",chr)
    if(is.null(main.title)) {
        main.title<-paste("Manhattan Plot for association between methylation and "
                ,x$info$Phenotype,sep="")
          }
 
    cpgtitle<-data.frame(cpglab=as.character(cpgname),chr,pos)
    
    info<-merge(test,cpgtitle,by.x="CPG.Labels",by.y="cpglab",sort=FALSE)
    score<-info$P.value
    cutoff<-length(score)-nrow(x$Holm.sig)+1
    if(cutoff>length(score)) cutoff=length(score)
    chr<-as.character(info$chr)
    pos<-info$pos
    chr[chr=="X"]=23
    chr[chr=="Y"]=24
    
    chr<-as.numeric(chr)
    pos<-as.numeric(pos)
    k=order(chr,pos)
    chr=chr[k]
    pos=pos[k]
    
    score=-log(score[k],base=10)
    cmax=max(chr)
    missingspots<-sum(is.na(chr))>0
    if(is.na(cmax)){
      probspots<-which(is.na(chr))
      oldmax<-max(chr,na.rm=T)
      chr[probspots]=oldmax+1
      cmax<-max(chr)
      topick<-range(pos[chr==oldmax])
      pos[probspots]<-seq(topick[1],topick[2],length.out=length(probspots))
      k=order(chr,pos)
      chr=chr[k]
      pos=pos[k]
    }

    genomepos=pos
    if(!is.null(chr.list)) {
        cpg.subset<-which(chr %in% chr.list)
        score<-score[cpg.subset]
        genomepos<-genomepos[cpg.subset]
        chr<-chr[cpg.subset]
           }
   chr_unique<-unique(chr)
    if(length(chr_unique)>1) {
    for (i in chr_unique[2:length(chr_unique)]) {
      genomepos[chr==i]=genomepos[chr==i]+max(genomepos[chr==chr_unique[which(chr_unique==i)-1]])+100
          }
    }

    
    bonsig<- which(score > -log(.05/cutoff,base=10))
    lessig<-which(score < -log(.05/cutoff,base=10) & score > -log(.05/cutoff,base=10)/2)
    chrcolor<-chr
  
  
    original<-par(no.readonly = TRUE)
    if(!is.null(color.list)){
      number.chrs<-length(chr_unique)
      unique.colors<-unique(chrcolor)
      if(number.chrs < length(color.list)){
          replacement.color<-color.list[1:number.chrs]
            }
      else {
        replacement.color<-c(rep(color.list,floor(number.chrs/length(color.list))))
        remain.cpg<-number.chrs%%length(color.list)
        if(remain.cpg!=0) {
          replacement.color<-c(replacement.color,color.list[1:remain.cpg])
            }
         }
        chrcolor.new=chrcolor
        for (i in 1:number.chrs) {
            chrcolor.new[which(chrcolor == unique.colors[i])] <- replacement.color[i]
        }
       chrcolor=chrcolor.new
      chrcolor<-as.character(chrcolor)
              }
     else{
        chrcolor[which(chrcolor=="7" | chrcolor=="15" | chrcolor=="23")]="orange"
      }
      y.lab<-c(0,max(-log(.05/cutoff,base=10),score,na.rm=TRUE)+.05)
    plot(range(genomepos),y.lab,type="n",xaxt='n',bty='7',xlab="chromosome",
          ylab=expression(paste("Observed -log ", scriptstyle(10), "(P-values)",sep="")),main=main.title,...)
    if(is.null(point.size)){
     points(genomepos,score,col=chrcolor,cex=.2)
     points(genomepos[bonsig],score[bonsig],col=chrcolor[bonsig],pch=16)
     points(genomepos[lessig],score[lessig],col=chrcolor[lessig],pch=16,cex=.5)
    }
    else{
     points(genomepos,score,col=chrcolor,cex=point.size,...)
     points(genomepos[bonsig],score[bonsig],col=chrcolor[bonsig],cex=point.size,...)
     points(genomepos[lessig],score[lessig],col=chrcolor[lessig],pch=16,cex=point.size,...)
    }

    

    abline(-log(.05/cutoff,base=10),0)
    if(nrow(x$FDR.sig)>0) {
       abline(-log(max(x$FDR.sig$P.value),base=10),0,lty=2)
               }
    if(!is.null(cpg.labels)) {
      
       if(cpg.labels=="FDR" ) {
        if(nrow(x$FDR.sig)>0){
          fdr.cutoff<- -log(max(x$FDR.sig$P.value),base=10)
          fdr.labels<-x$results$CPG.Labels[k][which(score>=fdr.cutoff)]
          text(genomepos[which(score>fdr.cutoff)],score[which(score>fdr.cutoff)],labels=fdr.labels,pos='4',cex=.7) 
          }
        else{
          warning("No labels displayed due to no FDR significant sites\n")
            }    
          }

       if(cpg.labels=="HOLM") {
        if(nrow(x$Holm.sig)>0){
          holm.labels<-x$results$CPG.Labels[k][which(score> -log(.05/cutoff,base=10))]
          text(genomepos[which(score> -log(.05/cutoff,base=10))], score[which(score> -log(.05/cutoff,base=10))],
              labels=holm.labels,pos='4',cex=.7)
                   }
        else{
          warning("No labels displayed due to no Holm significant sites\n")
            }           
                  }
            }
    meds=matrix(0,cmax)
    labs=matrix(NA,cmax)
    tix=matrix(0,cmax+1)
    labs_null=matrix(NA,cmax+1)
    for (i in 1:length(chr_unique)) {
      tix[i+1]=max(genomepos[chr==chr_unique[i]])
       meds[i]<-(tix[i]+tix[i+1,1])/2
      labs[i]=chr_unique[i] 
      }
    if(cmax>=23) {
      if(!missingspots){
        labs[which(labs==23)]="X"
        if(sum(chr==24)>0) {labs[which(labs==24)]="Y"} 
        }
      if(missingspots){
        if(cmax>23){
          labs[which(labs==23)]="X"
        }
        if(cmax>24){
          labs[which(labs==24)]="Y"
        }
        labs[which(labs==cmax)]="NoInfo"

      }
    }
    axis(side=1,at=meds,lwd.ticks=0,labels=labs)
    axis(side=1,at=tix,lwd.ticks=1,labels=labs_null)
        
    if(!is.null(save.plot))  {
      if(file.type=="eps") {
            dev.off()
            }
      else{
        if(popup.pdf) {
          dev.copy2pdf(file=paste(save.plot[1],".pdf",sep=""))
            }
        else{
            dev.off()
            }
      
      }}
       }
