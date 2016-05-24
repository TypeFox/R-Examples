#' Fire frequency using kernel density
#' 
#' Computes paleo-fire frequency for a set of fire events using a gaussian
#' kernel density estimation procedure based on a defined bandwidth (see Mudelsee 2004 for 
#' details). Pseudo-replicated values are used to correct for edge bias, equivalent to
#' "minimum slope" correction in Mann (2004).
#' 
#' @param fevent Numeric vector, set of dates
#' @param up Numeric, upper age for fire frequnecy calculus
#' @param lo Numeric, lower age for fire frequnecy calculus
#' @param bandwidth Numeric, bandwidth in years, or character for automatic
#' bandwidth calculation (e.g. "bw.ucv" for unbiased cross validation) see
#' \code{\link[stats]{bandwidth}} for details
#' @author O. Blarquez
#' @param nbboot Numeric, number of bootstrap replicates
#' @param alpha Numeric, confidence interval (default 0.01)
#' @return ff data.frame, with fire frequency, bandwidth and CIs
#' @seealso \code{\link{plot.kdffreq}}
#' @references  Mann, M. E. (2004). On smoothing potentially non-stationary
#' climate time series. Geophysical Research Letters, 31(7). \cr \cr 
#' Mudelsee, M., Börngen, M., Tetzlaff, G., & Grünewald, U. (2004). Extreme floods 
#' in central Europe over the past 500 years: Role of cyclone pathway “Zugstrasse Vb”.
#'  Journal of Geophysical Research: Atmospheres (1984–2012), 109(D23).
#' @examples
#' 
#' set.seed(123)
#'  fevent=c(round(abs(rnorm(20,mean=7,sd=5))*1000),round(abs(rnorm(10,mean=8,sd=1))*1000))
#' 
#'  ff=kdffreq(fevent,bandwidth = 1000, nbboot=10)

kdffreq=function(fevent,
                 up=NULL,
                 lo=NULL,
                 bandwidth=NULL,
                 nbboot=NULL,
                 alpha=NULL)
{    
  if(is.null(up)) up=min(fevent)
  if(is.null(lo)) lo=max(fevent)
  if(is.null(bandwidth)) bandwidth=1000
  if(is.null(alpha)) alpha=.01
  if(is.null(nbboot)) nbboot=999
  
  
  #Organize data
  fevent=fevent[!is.na(fevent)]
  fevent=sort(fevent);
  fevent_seq=lo-up;
  
  ## Pseudodata generation (only 20% of original data)
  
  # for the upper edge
  if (length(fevent[fevent>lo])!=0){
    pseudo_lo=fevent[fevent>lo]} else {
      pseudo_lo=-((fevent)-rep((max(fevent)),length(fevent)))+max(fevent)
      pseudo_lo=pseudo_lo[(length(pseudo_lo)-round(length(pseudo_lo)*0.2)):length(pseudo_lo)]
    }
  
  # for the lower edge
  if (length(fevent[fevent<up])!=0){
    pseudo_up=fevent[fevent<up] } else { #if data exist use it
      #elseif generate mirored events
      pseudo_up=-((fevent)-rep((min(fevent)),length(fevent)))+min(fevent)
      pseudo_up=pseudo_up[1:round(length(pseudo_up)*0.2)]
    }
  
  ## Bandwidth selection using "bandwidth" function (see stats and density)
  
  if (is.character(bandwidth))
  bandwidth=eval(parse(text=paste0(bandwidth,"(fevent)")))
  
  ## Kernel density estimation
  fevent_pseudo=sort(c(pseudo_up,fevent,pseudo_lo))
  n=length(fevent_pseudo); #number of events
  t=seq(min(fevent_pseudo), max(fevent_pseudo), 10); 
  
  # points where kernel estimates the
  # density
  ffreq1=density(fevent_pseudo,bw=bandwidth,kernel="gaussian",from=min(fevent_pseudo),
                 to=max(fevent_pseudo), n=length(t))
  # plot(ffreq1$x,ffreq1$y*n,type="l",ylim=c(0,0.01))
  
  ## Bootstrap
  fmat=matrix(nrow=length(t),ncol=nbboot)
  
  for(i in 1:nbboot){
  fmat[,i]=c(density(fevent_pseudo[unique(sample(1:n,replace=TRUE))],
                    bw=bandwidth,kernel="gaussian",from=min(fevent_pseudo),
                   to=max(fevent_pseudo), n=length(t))$y*n)
  #lines(t, fmat[,i],col="grey67")
  }
  
  CI=apply(fmat, 1, quantile, probs=c(alpha/2,1-(alpha/2)) )
  # lines(t,CI[1,])
  # lines(t,CI[2,])

  
  ff= data.frame(cbind(ffreq1$x,ffreq1$y*n,t(CI)))
  names(ff)=c("age","ff","lo","up")
  ff=ff[ff$age>up & ff$age<lo ,]
  class(ff)="kdffreq"
  return(ff)
}

#' plot.kdffreq
#' 
#' Plot fire frequency calculated using the \code{\link[paleofire]{kdffreq}} function 
#' 
#' @method plot kdffreq
#' @export
#' @param x Object returned by kdffreq
#' @param xlim Numeric x axis limits
#' @param ylim Numeric, y axis limits
#' @param ... other arguments
#' @seealso \code{\link[paleofire]{kdffreq}} 
#'  @examples
#'  
#'  set.seed(123)
#'  fevent=c(round(abs(rnorm(20,mean=7,sd=5))*1000),round(abs(rnorm(10,mean=8,sd=1))*1000))
#' 
#'  ff=kdffreq(fevent,bandwidth = 1000, nbboot=10)
#'  plot(ff)


plot.kdffreq=function(x,ylim=NULL,xlim=NULL,...){
  
  if(is.null(ylim)) ylim=c(min(x$lo),max(x$up))
  if(is.null(xlim)) xlim=c(min(x$age),max(x$age))
  
  plot(x$age,x$ff,type="l",xlab = "Age", ylab = "FF (#.yr-1)",
       xlim=xlim,ylim=ylim)
  lines(x$age,x$up)
  lines(x$age,x$lo)
  
  
}
