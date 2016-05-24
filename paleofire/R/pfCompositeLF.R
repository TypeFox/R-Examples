#' Produce a composite serie from multiple charcoal records using a local
#' regression procedure (from the locfit package)
#' 
#' Produces a composite series from multiple charcoal records by using a robust
#' locally weighted scatterplot smoother (LOWESS). The robust LOWESS uses the
#' locfit function from the locfit package and is applied repeatedly (nboot
#' times) on bootstrapped charcoal sites samples. The records charcoal values
#' are pre-binned prior to sites resampling. This procedure is equivalent to
#' Daniau et al. (2012).
#' 
#' 
#' @param TR An object returned by \code{\link[paleofire]{pfTransform}}
#' @param tarAge Numeric, the target ages for prebinning given in years (e.g.
#' tarAge = seq(0, 10000, 20)). If unspecified the sequence is defined as
#' tarAge=seq(from=min age, to=max Age, by=median resolution).
#' @param binhw Numeric, bin half width for the prebinning procedure (use the
#' same value as tarAge intervals for overlapping bins or tarAge intervals/2
#' for non-overlapping bins, defautl).
#' @param nboot Numeric, a number specifying the number of bootstrap
#' replicates.
#' @param hw Numeric, the half window width for the locfit procedure (in
#' years).
#' @param conf Numeric, define confidence levels.
#' @param pseudodata Logical, if TRUE 10 percent of the data is reflected at
#' the top and the bottom of the resampled serie prior of each locfit
#' regression in order to correct for the edge effect introduced by the local
#' regression, see Cowling & Hall (1996). Equivalent to "minimum slope"
#' correction in Mann(2004).
#' @param verbose Logical: verbose or not...
#' @return
#' 
#' \item{out}{A "pfCompositeLF" object.}
#' @author O.Blarquez
#' @references Daniau, A. L., P. J. Bartlein, S. P. Harrison, I. C. Prentice,
#' S. Brewer, P. Friedlingstein, T. I. Harrison-Prentice, J. Inoue, K. Izumi,
#' J. R. Marlon, S. Mooney, M. J. Power, J. Stevenson, W. Tinner, Andri, M., J.
#' Atanassova, H. Behling, M. Black, O. Blarquez, K. J. Brown, C. Carcaillet,
#' E. A. Colhoun, D. Colombaroli, B. A. S. Davis, D. D'Costa, J. Dodson, L.
#' Dupont, Z. Eshetu, D. G. Gavin, A. Genries, S. Haberle, D. J. Hallett, G.
#' Hope, S. P. Horn, T. G. Kassa, F. Katamura, L. M. Kennedy, P. Kershaw, S.
#' Krivonogov, C. Long, D. Magri, E. Marinova, G. M. McKenzie, P. I. Moreno, P.
#' Moss, F. H. Neumann, E. Norstrom, C. Paitre, D. Rius, N. Roberts, G. S.
#' Robinson, N. Sasaki, L. Scott, H. Takahara, V. Terwilliger, F. Thevenon, R.
#' Turner, V. G. Valsecchi, B. Vanniere, M. Walsh, N. Williams, and Y. Zhang.
#' 2012. Predictability of biomass burning in response to climate changes.
#' Global Biogeochem. Cycles 26:GB4007. \cr \cr Cowling A, Hall P (1996) On
#' pseudodata methods for removing boundary effects in kernel density
#' estimation. Journal of the Royal Statistical Society, Series B 58(3):
#' 551-563. \cr \cr Mann, M. E. (2004). On smoothing potentially non-stationary
#' climate time series. Geophysical Research Letters, 31(7).
#' @examples
#' 
#' ID=pfSiteSel(id_region=="WNA0", l12==1, long>=-160 & long<=-140)
#' plot(ID, xlim=c(-180, -130), ylim=c(40,80))
#' TR=pfTransform(ID, method=c("MinMax","Box-Cox","MinMax","Z-Score"),
#'                BasePeriod=c(200,2000),QuantType="INFL")
#' 
#' COMP1=pfCompositeLF(TR, tarAge=seq(-50,4000,10), hw=200, nboot=100)
#' 
#' plot(COMP1)
#' 
#' ## Note: comparing confidence intervals based on 100 replicates is not recommended
#' # (100 is used to decrease analysis time)
#' 
#' 
pfCompositeLF=function(TR,hw=250,
                       tarAge=NULL,binhw=NULL,
                       nboot=1000,conf=c(0.05,0.95),
                       pseudodata=FALSE,verbose=TRUE)
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
  
  ## Prebinning procedure
  # Define the sequence for binning if unspecified
  if(is.null(tarAge)){
    AgeRes=matrix(nrow=length(TR$Age[,1])-1,ncol=length(TR$Age[1,]))
    for (i in 1:length(TR$Age[1,])){
      AgeRes[,i]=c(diff(TR$Age[,i]))
    }
    width=ceiling(median(na.omit(AgeRes))/10)*10
    binI=floor(min(na.omit(TR$Age))/10)*10
    binF=ceiling(max(na.omit(TR$Age))/10)*10
    tarAge=seq(binI,binF,width)
    if(is.null(binhw))
      binhw=(tarAge[2]-tarAge[1])/2
  }
  
  m=length(TR$TransData[,1])
  n=length(TR$TransData[1,])
  
  # Matrix to store results
  result=matrix(ncol=length(TR$Age[1,]),nrow=length(tarAge))
  
  ## Use non-overlapping bins by default if unspecified
  if(is.null(binhw))
    binhw=(tarAge[2]-tarAge[1])/2
  
  ## Prebinning procedure

  if(verbose==TRUE){
    percent=seq(10,100,by=10)
    values=round(percent*n/100)
    cat("Prebinning...")
    cat("\n")
    cat("Percentage done: ")
  }
  
  for (k in 1:n){
    if(length(TR$TransData[is.na(TR$TransData[,k])==FALSE,k])!=0){
      for (i in 1:length(tarAge)){
        t=na.omit(cbind(as.numeric(TR$Age[,k]),as.numeric(TR$TransData[,k])))
        result[i,k]=mean(t[t[,1]>tarAge[i]-binhw & t[,1]<tarAge[i]+binhw,2])
      }
    }
    ## print
    if(k %in% values & verbose==TRUE)
      cat(percent[values==k]," ",sep="")
  }
  if(verbose==TRUE) 
    cat("\n")
  
  # suppres Inf values occuring with specific charcoal series (binary series)
  result[!is.finite(result)]=NA
  
  # Target Ages:
  centres=tarAge
  
  # Matrix to strore boot results
  mboot=matrix(nrow=length(centres),ncol=nboot)
  
  set.seed(1)
  
  if(verbose==TRUE){
    percent=seq(10,100,by=10)
    values=round(percent*nboot/100)
    cat("Bootstrapping...")
    cat("\n")
    cat("Percentage done: ")
  }
  
  ## Bootstrap procedure (with locfit)
  for (i in 1:nboot){
    ne=sample(seq(1,length(result[1,]),1),length(result[1,]),replace=TRUE)
    y=as.vector(result[,ne])
    x=as.vector(rep(centres,length(ne)))
    dat=na.omit(data.frame(x=x,y=y))
    #dat=dat[ order(dat$x), ]
    
    if (pseudodata==TRUE){
      ## Pseudodata part:
      # Generate a mirror of data series at top and bottom
      # Only reflect 10% of the data
      pseudo_n=round(length(dat$x)*0.3)
      
      ## Upper part
      pseudo_up=-((dat$x)-rep((min(dat$x)),length(dat$x)))+min(dat$x)
      pseudo_up=as.data.frame(cbind(pseudo_up,dat$y))
      pseudo_up=pseudo_up[1:pseudo_n,]
      pseudo_up=pseudo_up[length(pseudo_up[,1]):1,]
      colnames(pseudo_up)=c("x","y")
      ## Lower part
      pseudo_lo=-((dat$x)-rep((max(dat$x)),length(dat$x)))+max(dat$x)
      pseudo_lo=as.data.frame(cbind(pseudo_lo,dat$y))
      pseudo_lo=pseudo_lo[(length(dat$x)-pseudo_n):length(dat$x),]
      pseudo_lo=pseudo_lo[length(pseudo_lo[,1]):1,]
      colnames(pseudo_lo)=c("x","y")
      # New random serie with 10% pseudodata
      dat=rbind(pseudo_up,dat,pseudo_lo)
    }
    ##
    x=as.vector(dat$x)
    y=as.vector(dat$y)
    locboot <- locfit(y ~ lp(x, deg = 1, h = hw), maxk = 2000, family = "qrgauss")
    if(is.na( locboot$dp[7])==FALSE){
    predboot <- predict(locboot, newdata = centres, se.fit = TRUE)
    mboot[, i] <- predboot$fit}
    
    # Verbose
    if(i %in% values & verbose==TRUE)
      cat(percent[values==i]," ",sep="")
  }
  
  bootci=t(apply(mboot, 1, quantile, probs = conf,  na.rm = TRUE))
  bootmean=t(apply(mboot, 1, mean,  na.rm = TRUE))
  bootmed=t(apply(mboot, 1, median,  na.rm = TRUE))
  
  ## Locfit of all stacked data
  rm(x,y,dat)
  y=c(result)
  x=rep(centres,length(result[1,]))
  dat=as.data.frame(cbind(x,y))
  dat=na.omit(dat[ order(x), ])
  if (pseudodata==TRUE){
    ## Pseudodata part:
    # Generate a mirror of data series at top and bottom
    # Only reflect 10% of the data
    pseudo_n=round(length(dat$x)*0.3)
    
    ## Upper part
    pseudo_up=-((dat$x)-rep((min(dat$x)),length(dat$x)))+min(dat$x)
    pseudo_up=as.data.frame(cbind(pseudo_up,dat$y))
    pseudo_up=pseudo_up[1:pseudo_n,]
    pseudo_up=pseudo_up[length(pseudo_up[,1]):1,]
    colnames(pseudo_up)=c("x","y")
    ## Lower part
    pseudo_lo=-((dat$x)-rep((max(dat$x)),length(dat$x)))+max(dat$x)
    pseudo_lo=as.data.frame(cbind(pseudo_lo,dat$y))
    pseudo_lo=pseudo_lo[(length(dat$x)-pseudo_n):length(dat$x),]
    pseudo_lo=pseudo_lo[length(pseudo_lo[,1]):1,]
    colnames(pseudo_lo)=c("x","y")
    # New serie with 10% pseudodata
    dat=rbind(pseudo_up,dat,pseudo_lo)
  }
  x=as.vector(dat$x)
  y=as.vector(dat$y)
  locbootA <- locfit(y ~ lp(x, deg = 1, h = hw), maxk = 500, family = "qrgauss")
  predbootA <- predict(locbootA, newdata = centres, se.fit = TRUE)
  locfitAll <- predbootA$fit
  
  # write out the transformed data
  result2=as.data.frame(cbind(centres,locfitAll,t(bootmean),bootci))
  colnames(result2)=c("AGE","LocFit","MEAN(of_boot)",as.character(conf))
  
  colnames(result)=colnames(TR$TransData)
  
  colnames(result)=colnames(TR$TransData)
  output=structure(list(BinnedData=structure(result,row.names = as.character(centres),col.names=colnames(TR$TransData), class = "matrix"),
                        Result=result2,
                        mboot=mboot,
                        BinCentres=centres,    
                        BinWidth=binhw*2,
                        nboot=nboot,
                        halfwidth=hw,
                        conf=conf,
                        locfitAll=locfitAll,
                        BootMed=structure(bootmed,row.names = as.character(centres),col.names=c("Median"),class = "matrix"), 
                        BootMean=structure(bootmean,row.names = as.character(centres),col.names=c("Mean"),class = "matrix"), 
                        BootCi=structure(bootci,row.names = as.character(centres),class = "matrix") ))
  class(output)="pfCompositeLF"
  return(output)  
  if(verbose==TRUE) 
    cat("\n")
}

######SUMMARY########


######PLOT########





#' plot.pfCompositeLF
#' 
#' Plot pfCompositeLF object
#' 
#' @export
#' @method plot pfCompositeLF
#' @param x A "pfCompositeLF" object.
#' @param type Character, type of plot among "ci", "prctile", "density"
#' @param add Character, add=NULL by default, add="sitenum" could be specified
#' to plot the sites number in eah bin along with the composite curve.
#' @param conf Numeric, confidence levels.
#' @param palette Character, color palette used with type=c("prctile",
#' "density") among "jet" and "BW".
#' @param xlim Numeric, x axis limits.
#' @param ylim Numeric, y axis limits.
#' @param main Character, title of the plot.
#' @param text Logical, text options.
#' @param what Character, indicates which transformed charcoal trend is used
#' for the plot (type="ci"), default "locfit" indicates that the trend is the
#' locfit applied to All binned data, use "mean" or "median" to plot the mean
#' or median of the locfit relplicates given by the bootstrap procedure.
#' @param \dots \dots{}
#' @author O. Blarquez
#' @examples
#' 
#' ID=pfSiteSel(id_region=="WNA0",l12==1,long>=-160,long<=-140)
#' 
#' TR=pfTransform(ID, method=c("MinMax","Box-Cox","MinMax","Z-Score"),
#'                BasePeriod=c(200,2000),QuantType="INFL")
#' 
#' COMP1=pfCompositeLF(TR, tarAge=seq(-50,4000,10), hw=200, nboot=100)
#' 
#' plot(COMP1)
#' 
#' ## Note: comparing confidence intervals based on 100 replicates is not recommended
#' # (100 is used to decrease analysis time)
#' 
#' 
plot.pfCompositeLF=function(x,type="ci",add="NULL",conf=c(0.05,0.95),palette="jet",xlim=NULL,
                            ylim=NULL,main="Composite",text=FALSE,what="locfit",...){
  # Value for plotting:
  w=(x$BinCentres[2]-x$BinCentres[1])/2
  
  # Lims
  if (type=="ci") bootci1=t(apply(x$mboot, 1, quantile, probs = conf,  na.rm = TRUE)) 
  if (type=="prctile") bootci1=t(apply(x$mboot, 1, quantile, probs = seq(0, 1, .01),  na.rm = TRUE))
  if(is.null(xlim)) xlim=c(max(x$BinCentres)+w,min(x$BinCentres)-w)
  if(is.null(ylim)) {
    #bootci1=t(apply(x$mboot, 1, quantile, probs = seq(0, 1, .01),  na.rm = TRUE))
    #ylim= c(min(bootci1,na.rm=T),max(bootci1,na.rm=T))
    ylim= range( x$Result[,4:5])}
  
  # What to plot as the trend:
  if (what=="mean") trend=x$BootMean
  if (what=="median") trend=x$BootMed
  if (what=="locfit") trend=x$locfitAll
  
  # Now plot:
  
  if (type=="ci"){
    
    if(add=="sitenum")
      par(mfrow=c(2,1))
    
    bootci1=t(apply(x$mboot, 1, quantile, probs = conf,  na.rm = TRUE))    
    
    plot(x$BinCentres,trend, 
         xlim=xlim, ylim=ylim, axes=F, mgp=c(2,0,0),
         main=main, font.main=1, lab=c(8,5,5), 
         ylab="Composite", xlab="Age (cal yr BP)", cex.lab=1, pch=16, cex=0.5, type="l")
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
      plot(x$BinCentres,sitenum,xlim=c(max(x$BinCentres)+w,min(x$BinCentres)-w), 
           ylim= c(min(sitenum,na.rm=T),max(sitenum,na.rm=T)), axes=F, mgp=c(2,0,0),
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
    plot(NULL, type = "n", 
         xlim=xlim, 
         ylim=ylim,axes=FALSE,ylab="Composite",xlab="Age",main=main)
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
    #bootci1=t(apply(x$mboot, 1, quantile, probs = seq(0, 1, .01),  na.rm = TRUE))
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
    image(x$BinCentres, seqI, t(img),col = pal(100),xlab="Age",ylab="Composite",main=main,axes=F,
          xlim=xlim)
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

