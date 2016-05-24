plot.cpg.perm <-
function(x,save.plot=NULL,file.type="pdf",popup.pdf=FALSE,main.title=NULL,eps.size=c(5,5),tplot=FALSE,perm.ci=TRUE,classic=TRUE,
          gc.p.val=FALSE,gcdisplay=FALSE,...) {
  if(x$perm.p.values$nperm <  100 | !perm.ci) {
        plot.cpg(x,save.plot,tplot,...) }

  else {
  
      if(!is.null(save.plot)){
      if(!(file.type %in% c("eps","pdf"))) {
           stop("Incorrect file type. Must be pdf or eps\n")
               }
     
     
      if(file.type=="eps"){
        postscript(paste(save.plot[1],".eps",sep=""), horizontal = FALSE, 
                onefile = FALSE, paper = "special",width=eps.size[1],height=eps.size[2])
             }
      if(file.type=="pdf" & !popup.pdf) {
         pdf(file=paste(save.plot[1],".pdf",sep=""))
          }
        }
  ob=x$results[which(!is.na(x$results$P.value)),]
   if(!is.factor(x$indep)) {
      
      gcvalue<-format(median(ob$T.statistic**2)/0.4549364,digits = 3)
      }
      
 else{
     gcvalue<-format((nlevels(x$indep)-1)*median(ob$F.statistic)/qchisq(.5,nlevels(x$indep)-1),digits = 3)
        }
     
     
  gcvalue<-paste("Genomic control factor = ", gcvalue,sep="")
  
    if(gc.p.val) {
        ob$P.value<-ob$gc.p.value
        fdr.method<-x$info$FDR.method
        holm.adj<-p.adjust(ob$P.value,"holm")
     if (fdr.method=="qvalue") {
         if (!requireNamespace("qvalue", quietly = TRUE)) {
            stop("qvalue needed for this to work. Please install it.",
                      call. = FALSE)
                      }
        fdr.adj<-tryCatch(qvalue::qvalue(ob$P.value), error = function(e) NULL)
         if(is.null(fdr.adj)) {
          fdr.adj <- tryCatch(qvalue::qvalue(ob$P.value, pi0.method = "bootstrap"), error = function(e) NULL)
          if(is.null(fdr.adj)) {
            fdr.method="BH"
        }}

          }
    if(fdr.method!="qvalue") {
          fdr.adj<-p.adjust(ob$P.value,fdr.method)
          }
    x$FDR.sig<-subset(ob,fdr.adj<.05)
    x$Holm.sig<-subset(ob,holm.adj<.05)
    
    
          }
          
  
  
  u=(1:nrow(ob)-1/2)/nrow(ob)
  pvalues<-ob$P.value
  adjusp<-p.adjust(ob$P.value,"holm")
  nam.ind<-x$info$Phenotype
  sig<-nrow(x$FDR.sig)
    if(sig>0) {
      nonsig<-ob$P.value[-(which(ob$P.value %in% x$FDR.sig$P.value))]
      }
    else {
      nonsig<-ob$P.value
      }
  u2=(1:length(pvalues)-1/2)/length(pvalues)
  if(!classic) {
    u=(1:length(nonsig)-1/2)/length(nonsig)
    pvalues<-nonsig
        }
  if(!tplot | is.factor(x$indep)) {

     if(tplot & is.factor(x$indep)) {
      warning("Can not do t-statistic plot with a factor variable\n")
      tplot=FALSE  }
    if(is.null(main.title)) {
       main.title<-paste("QQ plot for association\nbetween methylation and",
            nam.ind)
          }
    if(gc.p.val) {
        main.title<-paste("GC adjusted P-values:",main.title)
          }
     ob<-ob$P.value
    k=order(ob)
    legend.text.use<-c("Holm-significant",paste("FDR-significant (",as.character(x$info$FDR.method),")"),"95% confidence interval",gcvalue)
    l.lty     <-c(-1,-1,2,NA)
    l.pch     <-c(19,1,-1,NA)
    l.col     <-c("red","red","black","black")
    if(!gcdisplay){
      legend.text.use<-legend.text.use[-4]
      l.lty          <-l.lty[-4]
      l.pch          <-l.pch[-4]
      l.col          <-l.col[-4]
    } 
    plot(-log(u,base=10),-log(sort(pvalues),base=10),xlab=expression(paste("Expected -log", scriptstyle(10), "(P-values)",sep="")),
            ylab=expression(paste("Observed -log ", scriptstyle(10), "(P-values)",sep="")),main=main.title,ylim=c(0,-log(min(ob),base=10)),cex=pointsizefunction(sort(pvalues)),...)
    
    
    legend(0,-log(min(ob),base=10),legend.text.use,lty=l.lty,pch=l.pch,col=l.col)
    
    if(sig>0) {
      if(nrow(x$Holm.sig) >0){
        holm<-sort(ob[which(adjusp<.05)])
        if(!classic) {
          points(rep(-log(min(u),base=10),length(holm)),-log(holm,base=10),col="red",pch=19)
            }
        if(classic) {
          points(-log(u[1:length(holm)],base=10),-log(holm,base=10),col="red",pch=19)
            }
        }
      if(nrow(x$Holm.sig)==0 | nrow(x$Holm.sig) < nrow(x$FDR.sig)) {
        fdrsig<-sort(ob[which(adjusp > .05 & ob %in% x$FDR.sig$P.value)])
        if(!classic) {
          points(rep(-log(min(u),base=10),length(fdrsig)),-log(fdrsig,base=10),col="red")
            }
        else {
           points(-log(u[(length(holm)+1):(length(holm)+length(fdrsig))],base=10),-log(fdrsig,base=10),col="red")
              }
      }}
    abline(a=0,b=1)
    lu=length(u)
    uu=1:lu
    quant<-x$perm.pval

    if(!classic & sig>0) {  
        quant<-quant[-((length(nonsig)+1):nrow(quant)),] 
      }

      points(-log(u,base=10),quant[,1],type='l',lty=2)
      points(-log(u,base=10),quant[,2],type='l',lty=2)

      }
  else {
    if(is.null(main.title)) {
       main.title<-paste("Expected vs. Observed T-statistics",
            "for association\nbetween methylation and",x$info$Phenotype)
          }
       if(gc.p.val) {
        main.title<-paste("GC adjusted T-statistics:",main.title)
          }
    
  
    k=order(ob$T.statistic)
    index3<-ob$T.statistic %in% x$FDR.sig$T.statistic

    index3<-index3[k]
    
    ob<-ob[k,]
    if(ncol(x$coefficients)==5) {
          df_use<-x$coefficients[!is.na(x$coefficients[,2]),2][k]
           }
    if(ncol(x$coefficients)==4) {
          df_use<-x$coefficients[!is.na(x$coefficients[,1]),1][k]
        }
    t.val<-qt(u2,df_use)

    adjusp<-adjusp[k]
    tstatistic<-ob$T.statistic[k]
    if(sig>0) { nonsig<-ob$T.statistic[-(which(ob$T.statistic %in% x$FDR.sig$T.statistic))]}
     else {nonsig<-ob$T.statistic}

    k2<-order(nonsig)
    remsig<-which(ob$T.statistic %in% x$FDR.sig$T.statistic)
  if(!classic) {
        tstatistic <-nonsig[k2]
        if(length(remsig)>0) {
             df_use<-df_use[-remsig][k2]
                }}
    y.val<-ob$T.statistic[k]
    t.val2<-qt(u,df_use)
    x.lim<-c(min(t.val,na.rm=TRUE),max(t.val,na.rm=TRUE))
    y.lim<-c(min(y.val,na.rm=TRUE),max(y.val,na.rm=TRUE))
    legend.text.use<-c("Holm-significant","Holm-significant",paste("FDR-significant (",as.character(x$info$FDR.method),")"),paste("FDR-significant (",as.character(x$info$FDR.method),")"),
                                    "95% confidence interval", gcvalue)
    l.lty     <-c(-1,-1,-1,-1,2,NA)
    l.pch     <-c(19,19,1,1,-1,NA)
    l.col     <-c("red","green","red","green","black","black")
    if(!gcdisplay){
      legend.text.use<-legend.text.use[-6]
      l.lty          <-l.lty[-6]
      l.pch          <-l.pch[-6]
      l.col          <-l.col[-6]
    } 
  
    qqplot(t.val2,tstatistic,xlab="Expected",ylab="Observed",main=main.title,xlim=x.lim,ylim=y.lim,...)


    legend(min(t.val,na.rm=TRUE),max(y.val,na.rm=TRUE),legend.text.use,lty=l.lty,pch=l.pch,col=l.col)
 
    if(sig>0) {
      holmreds<-which(adjusp<.05 & (y.val>0) )
      holmgreen<-which(adjusp<.05 & (y.val<0))
      if(nrow(x$Holm.sig) >0){
        if(!classic) {
          points(rep(max(t.val2),length(holmreds)),y.val[holmreds], col="red",pch=19)
          points(rep(min(t.val2),length(holmgreen)),y.val[holmgreen], col="green",pch=19)
           }
        else {
          if(length(holmreds)>0) {points(t.val2[(length(t.val2)-length(holmreds)+1):length(t.val2)],y.val[holmreds], col="red",pch=19)}
          if(length(holmgreen)>0) {points(t.val2[1:length(holmgreen)],y.val[holmgreen], col="green",pch=19)}
           }
        }
      if(nrow(x$Holm.sig)==0 | nrow(x$Holm.sig) < nrow(x$FDR.sig)) {
        sorted<-ob$P.value
        holmlered<-which(adjusp > .05 & (sorted %in% x$FDR.sig$P.value) & y.val >0)
        holmlegreen<-which(adjusp > .05 & sorted %in% x$FDR.sig$P.value & y.val<0)
        if(!classic){
          points(rep(max(t.val2),length(holmlered)),y.val[holmlered], col="red")
          points(rep(min(t.val2),length(holmlegreen)),y.val[holmlegreen], col="green")
          }
        else{
         if(length(holmlered)>0) { points(t.val2[(length(t.val2)-length(holmreds)-length(holmlered)+1):(length(t.val2)-length(holmreds))],y.val[holmlered],col="red") }
         if(length(holmlegreen)>0) { points(t.val2[(1+length(holmgreen)):(length(holmgreen)+length(holmlegreen))],y.val[holmlegreen],col="green")}
            }
      }
       }
    reds<-which((index3==TRUE)&(y.val>0))
    greens<-which((index3==TRUE)&(y.val<0))

      quant<-x$perm.tstat
     
      if(!classic & sig>0) {
        toremove<-c((nrow(quant)):(nrow(quant)+1-length(reds)),1:length(greens)) 
       quant<-quant[-toremove,]
       }
      points(t.val2,quant[,1],type='l',lty=2)
      points(t.val2,quant[,2],type='l',lty=2)

         }
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
      
      
      } }
