#' Produce a composite serie from multiple charcoal records
#' 
#' Produce a composite serie from multiple charcoal records using bootstrap
#' resampling, the sites charcoal values are binned and the mean in each bin is
#' calculated prior the bootstrap procedure. This procedure is equivalent to
#' Power et al. 2008.
#' 
#' 
#' @param TR An object returned by \code{\link{pfTransform}}
#' @param bins Numeric, the sequence for binning given in years (e.g.
#' bins=seq(from=0, to=10000, by=200)). If unspecified the sequence is defined
#' as bins=seq(from=min age, to=max age, by=median resolution).
#' @param nboot Numeric, a number specifying the number of bootstrap
#' replicates.
#' @param binning Logical, set to TRUE (default) for binning, if transformed
#' data are first interpolated this argument can be set to FALSE (no binning).
#' @param conf Numeric, define confidence levels.
#' @return Object of the class "pfComposite"
#' @author O.Blarquez
#' @references Power, M., J. Marlon, N. Ortiz, P. Bartlein, S. Harrison, F.
#' Mayle, A. Ballouche, R. Bradshaw, C. Carcaillet, C. Cordova, S. Mooney, P.
#' Moreno, I. Prentice, K. Thonicke, W. Tinner, C. Whitlock, Y. Zhang, Y. Zhao,
#' A. Ali, R. Anderson, R. Beer, H. Behling, C. Briles, K. Brown, A. Brunelle,
#' M. Bush, P. Camill, G. Chu, J. Clark, D. Colombaroli, S. Connor, A. L.
#' Daniau, M. Daniels, J. Dodson, E. Doughty, M. Edwards, W. Finsinger, D.
#' Foster, J. Frechette, M. J. Gaillard, D. Gavin, E. Gobet, S. Haberle, D.
#' Hallett, P. Higuera, G. Hope, S. Horn, J. Inoue, P. Kaltenrieder, L.
#' Kennedy, Z. Kong, C. Larsen, C. Long, J. Lynch, E. Lynch, M. McGlone, S.
#' Meeks, S. Mensing, G. Meyer, T. Minckley, J. Mohr, D. Nelson, J. New, R.
#' Newnham, R. Noti, W. Oswald, J. Pierce, P. Richard, C. Rowe, M. Sanchez
#' Goni, B. Shuman, H. Takahara, J. Toney, C. Turney, D. Urrego-Sanchez, C.
#' Umbanhowar, M. Vandergoes, B. Vanniere, E. Vescovi, M. Walsh, X. Wang, N.
#' Williams, J. Wilmshurst, and J. Zhang. 2008. Changes in fire regimes since
#' the Last Glacial Maximum: an assessment based on a global synthesis and
#' analysis of charcoal data. Climate Dynamics 30:887-907.
#' @examples
#' 
#' ## Composite charcoal record for North America:
#' ID=pfSiteSel(id_region=="WNA0" & l12==1)
#' plot(ID)
#' ## Transform data
#' res3=pfTransform(ID,method=c("MinMax","Box-Cox","Z-Score"),BasePeriod=c(200,4000))
#' 
#' ## Composite
#' comp=pfComposite(res3,bins=seq(from=0,to=12000,by=200))
#' plot(comp)
#' 
pfComposite=function(TR,
                     bins=NULL,
                     nboot=1000,
                     binning=TRUE,
                     conf=c(0.05,0.95))
{
  
  ## IF TR is a matrix
  if (is.matrix(TR) | is.data.frame(TR)){
    ID=unique(TR[,1])
    lgth=c()
    for (i in 1:length(ID)){
      lgth[i]=length(na.omit(TR[TR[,1]==ID[i],1]))
    }
    m=max(lgth)
    Age=matrix(nrow=m,ncol=length(ID))
    TransData=matrix(nrow=m,ncol=length(ID))
    for (i in 1:length(ID)) {
      Age[,i]=c(TR[TR[,1]==ID[i],3], rep(NA, m-length(TR[TR[,1]==ID[i],3])))
      TransData[,i]=c(TR[TR[,1]==ID[i],4], rep(NA, m-length(TR[TR[,1]==ID[i],4])))    
    } 
    colnames(TransData)=ID
    colnames(Age)=ID
    TR=structure(list(Age=structure(Age,class="matrix") ,TransData=structure
                      (TransData,class="matrix"),Method="unspecified"))
    
  }
  
  if(binning==TRUE){
    # Define the sequence for binning if unspecified
    if(is.null(bins)){
      AgeRes=matrix(nrow=length(TR$Age[,1])-1,ncol=length(TR$Age[1,]))
      for (i in 1:length(TR$Age[1,])){
        AgeRes[,i]=c(diff(TR$Age[,i]))
      }
      width=ceiling(median(na.omit(AgeRes))/10)*10
      binI=floor(min(na.omit(TR$Age))/10)*10
      binF=ceiling(max(na.omit(TR$Age))/10)*10
      bins=seq(binI,binF,width)
    }
    
    # If specified, values for plotting are:
    if(is.null(bins)==FALSE){
      width=bins[2]-bins[1]
      binI=bins[1]
      binF=bins[length(bins)]
    }
    
    # Matrix to store results
    result=matrix(ncol=length(TR$Age[1,]),nrow=length(bins)-1)
    #sm_result=matrix(ncol=length(ID),nrow=length(bins)-1)
    
    ## Binning procedure
    for (k in 1:length(TR$Age[1,])){
      c1 <- cut(TR$Age[,k], breaks = bins)
      tmean=tapply(TR$TransData[,k], c1, mean)
      result[,k]=c(as.numeric(tmean))
    } 
    # suppres Inf values occuring with specific charcoal series (binary series)
    result[!is.finite(result)]=NA
  }
  
  ##
  if(binning==FALSE){
    bins=TR$Age[,1]
    binF=max(bins)
    binI=min(bins)
    width=bins[2]-bins[1]
    mboot=matrix(nrow=length(bins),ncol=nboot)
    result=TR$TransData
  } else mboot=matrix(nrow=length(bins)-1,ncol=nboot)
  
  ## Bootstrap procedure
  for (i in 1:nboot){
    ne=sample(seq(1,length(result[1,]),1),length(result[1,]),replace=TRUE)
    ## If has one bin
    if(dim(result)[1]==1){
      mboot[1,i]=mean(result[,ne],na.rm=TRUE)
    } else mboot[,i]=c(apply(result[,ne],1,mean,na.rm=TRUE))
  }
  bootci=t(apply(mboot, 1, quantile, probs = conf,  na.rm = TRUE))
  bootmean=t(apply(mboot, 1, mean,  na.rm = TRUE))
  
  # write out the transformed data
  
  if(binning==TRUE){centres=bins[1:length(bins)-1]+width/2} else centres=bins
  if(binning==FALSE){width=NA}
  
  result2=as.data.frame(cbind(centres,t(bootmean),bootci))
  colnames(result2)=c("AGE","MEAN",as.character(conf))
  colnames(result)=colnames(TR$TransData)
  
  output=structure(list(BinnedData=structure(result,row.names = as.character(centres),col.names=colnames(TR$TransData), class = "matrix"),
                        mboot=mboot,
                        BinCentres=centres,    
                        Result=result2,
                        BinWidth=width,
                        nboot=nboot,
                        binning=binning,
                        BootMean=structure(bootmean,row.names = as.character(centres),col.names=c("Mean"),class = "matrix"), 
                        BootCi=structure(bootci,row.names = as.character(centres),col.names=as.character(conf),class = "matrix") ))
  class(output)="pfComposite"
  return(output)  
}

######SUMMARY########

######PLOT########





#' plot.pfComposite
#' 
#' Plot a pfComposite object
#' 
#' @method plot pfComposite
#' @export
#' @param x A "pfComposite" object.
#' @param type Character, type of plot among "ci", "prctile", "density"
#' @param conf Numeric, confidence levels.
#' @param palette Character, color palette used with type=c("prctile",
#' "density") among "jet" and "BW".
#' @param add Character, add="NONE" by default, add="sitenum" could be
#' specified to plot the sites number in eah bin along with the composite
#' curve.
#' @param main Character, title of the plot.
#' @param text Logical, text options.
#' @param \dots \dots{}
#' @author O. Blarquez
#' @examples
#' 
#' ## Composite charcoal record for North America:
#' ID=pfSiteSel(id_region=="WNA0",l12=1)
#' ## Transform data
#' res3=pfTransform(ID,method=c("MinMax","Box-Cox","Z-Score"),BasePeriod=c(200,4000))
#' 
#' ## Composite
#' comp=pfComposite(res3,bins=seq(0,5000,200))
#' plot(comp,type="density",smoothing=TRUE,spar=0.3)
#' 
plot.pfComposite=function(x,type="ci",conf=c(0.05,0.95),palette="jet",add="NONE",text=FALSE,main=NULL,...){
  # Value for plotting:
  w=(x$BinCentres[2]-x$BinCentres[1])/2
  
  if (type=="ci"){
    
    if(add=="sitenum")
    par(mfrow=c(2,1))
   
    bootci1=t(apply(x$mboot, 1, quantile, probs = conf,  na.rm = TRUE))    
    
    plot(x$BinCentres,x$BootMean, xlim=c(max(x$BinCentres)+w,min(x$BinCentres)-w), ylim= c(min(bootci1,na.rm=T),max(bootci1,na.rm=T)), axes=F, mgp=c(2,0,0),
         main=main, font.main=1, lab=c(8,5,5), 
         ylab="Composite", xlab="Age (cal yr BP)", cex.lab=1, pch=16, cex=0.5, type="o")
    axis(1); axis(2, cex.axis=1)
    axis(side = 1, at = seq(0, 99000, by = 500), 
         labels = FALSE, tcl = -0.2)  
    for (i in 1:length(conf)){
      lines(x$BinCentres,bootci1[,i],lty=2)
      pos=which.min(is.na(bootci1[,i]))
      if(text==TRUE) text(min(x$BinCentres)-200,bootci1[pos,i],paste(conf[i]*100,"%",sep=""),col="black")
    }
    
    # Plot site number
    if(add=="sitenum"){
    sitenum=length(x$BinnedData[1,])-rowSums(is.na(x$BinnedData))
    plot(x$BinCentres,sitenum,xlim=c(max(x$BinCentres)+w,min(x$BinCentres)-w), ylim= c(min(sitenum,na.rm=T),max(sitenum,na.rm=T)), axes=F, mgp=c(2,0,0),
         main=paste("Sites #"), font.main=1, lab=c(8,5,5), 
         ylab="Sites #", xlab="Age", cex.lab=0.8, pch=16, cex=0.5, type="o")
    axis(1); axis(2, cex.axis=1)
    axis(side = 1, at = seq(0, 99000, by = 500), 
         labels = FALSE, tcl = -0.2)
    }
  }
  
  
  
  if (type=="prctile"){
    bootci1=t(apply(x$mboot, 1, quantile, probs = seq(0, 1, .01),  na.rm = TRUE))
    bins1=x$BinCentres[is.na(bootci1[,1])==FALSE]
    bootci1=bootci1[is.na(bootci1[,1])==FALSE,]
    n=length(bootci1[1,])
    ## PLOT
    plot(NULL, type = "n", xlim=c(max(x$BinCentres)+w,min(x$BinCentres)-w), ylim = range(bootci1),axes=FALSE,ylab="Composite",xlab="Age",main=main)
    if(palette=="jet"){pal = colorRampPalette(rev(c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000")))} 
    if(palette=="BW"){
      pal = colorRampPalette(rev(c("white","black")))}
    xx=cbind(bins1,rev(bins1))
    coli=pal(50)
    for (i in 1:floor(n/2)) {
      yy <- cbind(as.vector(bootci1[,i]), rev(as.vector(bootci1[, n - i + 1])))
      polygon(xx, yy, col =coli[floor(n/2) - i + 1], border =coli[floor(n/2) - i + 1])  
    }
    for (i in c(2,11,51,91,100)){
      lines(bins1,bootci1[,i],col="grey",lty=2)
      if(text==TRUE) text(min(bins1)-200,median(bootci1[1:round(length(bootci1[,1])*0.02),i]),paste(i-1,"%",sep=""),col="grey")
    }
    axis(1)
    axis(side = 1, at = seq(0, 99000, by = 500), 
         labels = FALSE, tcl = -0.2) 
    axis(2)
  }
  
  if (type=="density"){
    seqI=seq(min(na.omit(x$mboot)),max(na.omit(x$mboot)),len=1000)
    img=matrix(nrow=1000,ncol=length(x$mboot[,1]))
    
    id=seq(1,length(x$mboot[,1]),1)[rowSums(x$mboot,na.rm=T)!=0]
    id[rowSums(x$mboot,na.rm=T)!=0]
    
    for (i in seq(1,length(x$mboot[,1]),1)[rowSums(x$mboot,na.rm=T)!=0]){
      kd=density(x$mboot[i,],na.rm=T)
      img[,i]=c(approx(kd$x,kd$y,seqI)$y)
    }
    layout(matrix(c(1,1,1,2), 1, 4, byrow = TRUE))
    if(palette=="jet"){pal = colorRampPalette((c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000")))} 
    if(palette=="BW"){
      pal = colorRampPalette((c("white","black")))}
    image(x$BinCentres, seqI, t(img),col = pal(100),xlab="Age",ylab="Composite",main=main,axes=F, xlim=c(max(x$BinCentres)+w,min(x$BinCentres)-w))
    axis(1, cex.axis=1, xaxp=c(0,99000,99)); axis(2, cex.axis=1)
    lines(x$BinCentres,rowMeans(x$mboot, na.rm=T))
    z=matrix(1:100,nrow=1)
    x=1
    y=seq(min(img,na.rm=T),max(img,na.rm=T),len=100)
    image(x,y,z,col=pal(100),axes=FALSE,xlab="Density",ylab="")
    axis(2)
    axis(side = 1, at = seq(0, 99000, by = 500), 
         labels = FALSE, tcl = -0.2) 
  }
}
