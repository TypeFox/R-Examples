#' Are cores sampled contiguously?
#' 
#' The function checks wether cores have been sampled contiguously or with a
#' depth resolution <1cm.
#' 
#' 
#' @param x An object of the class "pfSiteSel"
#' @param threshold Numeric, threshold for considering two samples as
#' contiguous (default=1cm)
#' @return Summary table of sites with the added contiguous logical column
#' (TRUE--FALSE)
#' @author O. Blarquez
#' @seealso \code{\link[paleofire]{pfResolution}}
#' @examples
#' 
#' x=pfSiteSel(lat>12,lat<60,long<(-50),long>-140)
#' contiguous(x)
#' 
contiguous=function(x,threshold=1){
  
  #x=pfSiteSel(lat>12,lat<60,long<(-50),long>-140)
  #plot(x)
  
  dat=pfExtract(x)
  
  # Keep preferably INFL units
  for(i in x$id_site){ 
    if(length(unique(dat[dat$ID_SITE==i,"TYPE"]))>1){
      un=unique(dat[dat$ID_SITE==i,"TYPE"])
      if("INFL" %in% un){
        dat[dat$ID_SITE==i  & dat$TYPE!="INFL","TYPE"]=NA
      } else {dat[dat$ID_SITE==i  & dat$TYPE!=un[1],"TYPE"]=NA}
    }
  }
  
  # x$id_site
  dat=na.omit(dat)
  
  # Contiguous???
  z=c()
  for(i in 1:length(x$id_site)){
    ddiff=diff(dat[dat$ID_SITE==x$id_site[i],2]*100)
    if(sum(ddiff,na.rm=T)!=0){
      z[i]=max(quantile(ddiff,probs=c(0.5,0.95)))} else {z[i]=0}
  }
  
  zz=z
  zz[zz>threshold+0.1]=0
  zz=as.logical(zz)
  sumsit=cbind(summary(x),contiguous=zz)
  class(sumsit)<-"contiguous"
  return(sumsit)
}





#' Plot "contiguous" object
#' 
#' Plot an object returned by contiguous, plot contiguous cores (or sites) in
#' green (T) and non-contiguous cores in red (F).
#' 
#' @method plot contiguous
#'
#' @param x An object returned by contiguous
#' @param ylim Numeric, ylim for the graph
#' @param xlim Numeric, xlim for the graph
#' @param \dots \dots{}
#' @return A plot.
#' @author O. Blarquez
#' @seealso \code{\link{contiguous}}
#' @examples
#' 
#' x=pfSiteSel(lat>12,lat<60,long<(-50),long>-140)
#' cont=contiguous(x)
#' plot(cont)
#' 
plot.contiguous=function(x,ylim=NULL,xlim=NULL,...){
  
  # sumsit
  coast=NULL; rm(coast)
  data(coast,envir = environment())
  
  if(is.null(xlim)) xlim=range(x$long)
  if(is.null(ylim)) ylim=range(x$lat)
  
  plot(coast$X,coast$Y,type='l',xlim=xlim,ylim=ylim)
  points(x$long[x$contiguous==T],x$lat[x$contiguous==T],cex=3,pch=16,col=rgb(0,1,0,0.5))
  points(x$long[x$contiguous==F],x$lat[x$contiguous==F],cex=3,pch=16,col=rgb(1,0,0,0.5))
  
  text(x$long[x$contiguous==T],x$lat[x$contiguous==T],"T")
  text(x$long[x$contiguous==F],x$lat[x$contiguous==F],"F")
}
