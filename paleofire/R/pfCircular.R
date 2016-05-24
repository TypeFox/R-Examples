#' Circular block bootstrap procedure applied to charcoal records compositing
#' results
#' 
#' Block bootstrap has been proposed to test the significances of changes in
#' stationary time series (Kunsch 1989). This procedure consists of splitting
#' each charcoal series into n-b+1 overlapping blocks of data, where n is
#' sample size and b the block size. These blocks are used to reconstruct
#' resampled individual charcoal series that are in turn used to estimate the
#' confidence intervals around the charcoal series composite mean.
#' 
#' 
#' @param comp A "pfComposite" object
#' @param b A numeric giving block size, if NULL the optimal block size for a
#' given series is given by: b= 2x(-1 /log(p)), where p is the lag one
#' autocorrelation coefficient of that series (Adams, Mann & Ammann 2003).
#' @param conf Numeric, calculated confidence intervals.
#' @param nboot Numeric, number of bootstrap replicates.
#' @param AgeLim Numeric, years defining a period to restrict the analysis to.
#' @return
#' 
#' \item{out}{A "pfCircular" object with estimated confidence intervals.}
#' @author O. Blarquez
#' @references Kunsch, H. R. 1989. The jackknife and the bootstrap for general
#' stationary observation s. The Annals of Statistics 17:1217-1241.
#' 
#' Adams, J. B., M. E. Mann, and C. M. Ammann. 2003. Proxy evidence for an El
#' Nino-like response to volcanic forcing. Nature 426:274-278.
#' @examples
#' 
#' ID=pfSiteSel(lat>49, lat<75, long>6, long<50)
#' plot(ID,zoom="world")
#' TR1=pfTransform(ID, method=c("MinMax","Box-Cox","Z-Score"),BasePeriod=c(200,2000))
#' 
#' ## Circular block bootstrapp
#' 
#' COMP=pfComposite(TR1, binning=TRUE, bins=seq(0,2000,100))
#' circ=pfCircular(COMP,conf=c(0.005,0.025,0.975,0.995),nboot=100)
#' plot(circ)
#' 
pfCircular=function(comp,b=NULL,conf=c(0.05,0.95),nboot=1000,AgeLim=NULL){
  
  ## R function developped from SEA.m   
  set.seed(1)

  ## Load matrix
  Temp=comp$BinnedData
  
  ## Define Age limits
  if(is.null(AgeLim)==FALSE){
    Temp=Temp[comp$BinCentres>AgeLim[1] & comp$BinCentres<AgeLim[2], ]
  }
  
  ## Block size calculus
  if(is.null(b)==TRUE){
    b=c()
    for(i in 1:length(Temp[1,])){
      r=cor(Temp[1:length(Temp[,1])-1,i],Temp[2:length(Temp[,1]),i],use="pairwise.complete.obs")
      yb=2*(-1/log(abs(r)))
      b[i]=c(ceiling(yb/mean(diff(comp$BinCentres))))
    }
    b[b==0 | b==1 | is.na(b) | is.finite(b)==FALSE]=2} else {b=rep(b,length(Temp[1,]))}
  
  
  ## Arrange data
  a=matrix(nrow=max(na.omit(b)),ncol=ncol(Temp))

  a[is.na(a)]=-999
  T_=rbind(a,Temp,a)
  
  ## Declare values for the boot process
  y_m=matrix(ncol=nboot,nrow=nrow(Temp))
  cat("# of Bootstrap:")
  
  for(k in 1:nboot){
    y_n=matrix(nrow=nrow(Temp),ncol=ncol(Temp))
    for(i in 1:length(Temp[1,])){
      n=ceiling(length(T_[,1])/b[i])
      q=trunc(length(T_[,1])-b[i])+1
      y=matrix(nrow=b[i],ncol=q)
      for(j in 1:q){
        y[,j]=c(T_[seq((j-1)+1,(j-1)+b[i],1),i])
      }
      o=sample(1:q,n*2,replace=TRUE)
      yy=c(y[,o])
      yy=yy[!yy == -999] 
      y_n[,i]=yy[1:length(Temp[,1])]
    }
    y_m[,k]=rowMeans(y_n,na.rm=TRUE)
    
    if(k %in% seq(0,nboot,10)) cat("", k)
  }
  
  ## Compile conf intervals
  boots=t(apply(y_m, 1, quantile, probs = conf,  na.rm = TRUE))
  
  ## Values for output
  if(is.null(AgeLim)==FALSE){
    yr=comp$BinCentres[comp$BinCentres>AgeLim[1] & comp$BinCentres<AgeLim[2]]
    Ci=comp$BootCi[comp$BinCentres>AgeLim[1] & comp$BinCentres<AgeLim[2], ]
    BootMean=comp$BootMean[comp$BinCentres>AgeLim[1] & comp$BinCentres<AgeLim[2]]
  } else {
    yr=comp$BinCentres
    Ci=comp$BootCi
    BootMean=comp$BootMean
  }
  
  ## Output
  output=structure(list(BootCirc=structure(boots,row.names = as.character(yr),col.names=conf, class = "matrix"),
                        conf=conf,
                        yr=yr,
                        BootCi=Ci,
                        BootMean=BootMean))
  class(output)="pfCircular"
  return(output)  
}

##--------------------------------------------------------------------------------------------------------




#' plot.pfCircular
#' 
#' Plot circular block bootstrap percentiles.
#' 
#' @method plot pfCircular
#' @export
#' @param x A "pfCircular" object.
#' @param ylim Numeric, x axis limits.
#' @param xlim Numeric, y axis limits.
#' @param ylab Character, y axis label.
#' @param xlab Character, x axis label.
#' @param main Character, title of the plot.
#' @param text Logical, text options.
#' @param \dots \dots{}
#' @author O. Blarquez
#' @examples
#' 
#' ID=pfSiteSel(lat>49,lat<75,long>6,long<50)
#' TR1=pfTransform(ID, method=c("MinMax","Box-Cox","Z-Score"),BasePeriod=c(200,2000))
#' 
#' ## Circular block bootstrapp
#' COMP=pfComposite(TR1, binning=TRUE, bins=seq(0,2000,100))
#' circ=pfCircular(COMP,conf=c(0.005,0.025,0.975,0.995),nboot=100)
#' plot(circ)
#' 
plot.pfCircular=function(x,ylim=NULL,xlim=NULL,ylab=NULL,xlab=NULL,main=NULL,text=FALSE,...){
  ## Plot
  
  t=c(x$BootMean,x$BootCirc)
  if(is.null(ylim)) ylim=c(min(t,na.rm=TRUE),max(t,na.rm=TRUE))
  if(is.null(xlim)) xlim=c(max(x$yr),min(x$yr))
  if(is.null(xlab)) xlab="Age (cal yr BP)"
  if(is.null(ylab)) ylab="Composite"
  if(is.null(main)) main=""
  
  plot(x$yr,x$BootMean,type="o",
       ylim=ylim,
       xlim=xlim,
       xlab=xlab,
       ylab=ylab, lab=c(8,5,5), 
       pch=16, cex=0.5,axes=F, mgp=c(2,0,0), main=main)
  for (i in 1:length(x$BootCirc[1,])){
    lines(x$yr,x$BootCirc[,i],lty=2)
    if(text==TRUE) text(min(x$yr)-200,x$BootCirc[1,i],paste(x$conf[i]*100,"%",sep=""),col="black")
  }
  axis(1); axis(2, cex.axis=1)
  axis(side = 1, at = seq(0, 99000, by = 500), 
       labels = FALSE, tcl = -0.2)  
}




