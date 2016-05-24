#' Power analysis plot of chi-squared test
#' 
#' @param es              effect size.
#' @param power           power of study
#' @param df              degree of freedom
#' @param sig.level       significance level
#' @param allele          in genetic association study, whether test allele or genotype
#' @param xlab            a title for the x axis
#' @param ylab            a title for the y axis
#' @param main            an overall title for the plot
#' @param grid            add grid lines or not
#' @param type            "np": plot sample size vs. power; "ne": plot effevct size vs. sample size
#' @seealso               \code{\link{power.chisq}}
#' @export
#' @examples
#' ## 'ne' type
#' ### multiple effect size and multiple power
#' es=seq(from=0.1,to=0.5,by=0.1);
#' power=seq(from=0.7,to=0.9,by=0.1);
#' power.plot.chisq(es=es,power=power,df=1,sig.level=0.05,type="ne")
#' power.plot.chisq(es=es,power=power,df=1,sig.level=0.05,type="np")
#' 
#' ### multiple effect size and single power
#' power.plot.chisq(es=seq(0.05,0.3,0.05),power=0.8,df=1,sig.level=0.05,type="ne")
#' power.plot.chisq(es=seq(0.05,0.3,0.05),power=0.8,df=1,sig.level=0.05,type="np")
#' 
#' ### single effect size and single power
#' power.plot.chisq(es=0.2,power=0.8,df=1,sig.level=0.05,type="ne")
#' power.plot.chisq(es=0.2,power=0.8,df=1,sig.level=0.05,type="np")
#' 
#' ### single effect size and multiple power
#' power.plot.chisq(es=0.2,power=seq(0.5,0.9,0.1),df=1,sig.level=0.05,type="ne")
#' power.plot.chisq(es=0.2,power=seq(0.5,0.9,0.1),df=1,sig.level=0.05,type="np")
power.plot.chisq <- function(es=NULL,power=NULL,df=NULL,sig.level=NULL,allele=FALSE,xlab=NULL, ylab=NULL,main=NULL,grid=FALSE,type=c("np","ne")){
  rr=NULL
  nes=length(es)
  npower=length(power)
  type=match.arg(type)
  
  n <- matrix(rep(0,nes*npower), nrow=nes,byrow=T)
  for (i in 1:npower){
    for (j in 1:nes){
      result <- power.chisq(es=es[j],df=df,n=NULL,power=power[i],sig.level=sig.level)
      if(allele==TRUE){
        n[j,i] <- ceiling(result$n/2)
      }else{
        n[j,i] <- ceiling(result$n)
      }
    }
  }
  
  if(type == "ne"){
    if(is.null(xlab)){xlab="Effect Size"}
    if(is.null(ylab)){ylab="Sample Size"}
    if(is.null(main)){main="Sample Size Estimation"}
    
    xrange <- range(es)
    yrange <- round(range(n))
    colors <- rainbow(length(power))
    
    plot=plot(xrange, yrange, type="n",xaxt="n",xlab=xlab,ylab=ylab,main=main)
    
    for (i in 1:npower){
      if(nes>1){
        lines(es, n[,i], type="l", lwd=2, col=colors[i])
      }else{
        points(es, n[,i], pch=16, col=colors[i])
      }
    }
    
    if(grid==TRUE){
      abline(v=0, h=seq(0,yrange[2],length.out=21), lty=2, col="grey89")
      abline(h=0, v=seq(xrange[1],xrange[2],0.01), lty=2,col="grey89")
    }
    
    axis(side=1,at=es,labels=es)
    legend("topright", title="power", as.character(power), fill=colors)
  }else if (type=="np"){
    if(is.null(xlab)){xlab="Sample size"}
    if(is.null(ylab)){ylab="Power"}
    if(is.null(main)){main="Sample Size Estimation"}
    
    xrange <- round(range(n))
    yrange <- range(power)
    colors <- rainbow(length(es))
    
    plot=plot(xrange, yrange, type="n",xlab=xlab,ylab=ylab,main=main)
    
    for (i in 1:nes){
      if(npower>1){
        lines(n[i,], power, type="l", lwd=2, col=colors[i])
      }else{
        points(n[i,],power, pch=16, col=colors[i])
      }
    }
    
    if(grid==TRUE){
      abline(v=0, h=seq(0,yrange[2],length.out=21), lty=2, col="grey89")
      abline(h=0, v=seq(xrange[1],xrange[2],0.01), lty=2,col="grey89")
    }
    
    #axis(side=1,at=n,labels=n)
    legend("bottomright", title="Effect Size", as.character(es), fill=colors, cex=0.6)
  }
  
  rownames=paste("effect.size=",es,sep="")
  colnames=paste("power=",power,sep="")
  dimnames(n)=list(rownames,colnames)
  rr$sig.level=sig.level
  rr$df=df
  rr$sample.size=n
  return(rr)
}

