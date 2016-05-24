pearsonDiagram <- function(max.skewness=sqrt(14),max.kurtosis=24,
                           squared.skewness=TRUE,
                           lwd=2,legend=TRUE,n=301) {
  x <- NaN                                                                      # make WARNING go away...
  t3fn  <- function(x) 3+3/2*x
  t5fn  <- function(x) -3*(16+13*x+2*sqrt(64+48*x+12*x^2+x^3))/(x-32)
  tnofn <- function(x) 1+x
  if (squared.skewness) {
    plot(0,0,type="n",xlim=c(0,max.skewness^2),ylim=c(0,max.kurtosis),xaxs="i",
         yaxs="i",xlab=expression(beta[1]),ylab=expression(beta[2]))
    b1seq <- seq(0,max.skewness^2,length.out=n)

    tmp <- list(x=c(0,max.skewness^2),y=tnofn(c(0,max.skewness^2)))             # type I
    if ((max(tmp$y)<max.kurtosis)&&(t3fn(max.skewness^2)>max.kurtosis)) {
      tmp$x <- c(tmp$x,max.skewness^2)
      tmp$y <- c(tmp$y,max.kurtosis)
    }
    tmp$x <- c(tmp$x,c(max.skewness^2,0))
    tmp$y <- c(tmp$y,t3fn(c(max.skewness^2,0)))
    polygon(tmp,col="yellow")

    tmp <- list(x = b1seq[b1seq<32], y = t5fn(b1seq[b1seq<32]))                 # type IV
    if (max(tmp$y)<max.kurtosis) {
      tmp$x <- c(tmp$x,max.skewness^2)
      tmp$y <- c(tmp$y,max.kurtosis)
    }
    tmp$x <- c(tmp$x,0); tmp$y <- c(tmp$y,max.kurtosis)
    polygon(tmp,col="lightcyan")

    tmp <- list(x=c(0,max.skewness^2),y=t3fn(c(0,max.skewness^2)))              # type VI
    if ((max(tmp$y)<max.kurtosis)&&(t5fn(max.skewness^2)>max.kurtosis)) {
      tmp$x <- c(tmp$x,max.skewness^2)
      tmp$y <- c(tmp$y,max.kurtosis)
    }
    tmp$x <- c(tmp$x, rev(b1seq[b1seq<32]))
    tmp$y <- c(tmp$y, t5fn(rev(b1seq[b1seq<32])))
    polygon(tmp,col="pink")

    curve(t5fn(x),from=0,to=min(31.999,max.skewness^2),n=n,add=TRUE,col="blue", # type V
          lwd=lwd)     
    lines(c(0,max.skewness^2),c(1,1+max.skewness^2),col="black",lwd=lwd)        # no Pearson
    lines(c(0,0),c(1,3),col="orange",lwd=lwd,xpd=NA)                            # type II
    lines(c(0,0),c(3,max.kurtosis),col="cyan",lwd=lwd,xpd=NA)                   # type VII
    lines(c(0,max.skewness^2),c(3,3+3/2*max.skewness^2),col="red",lwd=lwd)      # type III
    points(0,3,pch=16,col="brown",xpd=NA)                                       # type 0
  } else {
    plot(0,0,type="n",xlim=c(0,max.skewness),ylim=c(0,max.kurtosis),
         xaxs="i",yaxs="i",xlab=expression(sqrt(beta[1])),
         ylab=expression(beta[2]))
    b1seq <- seq(0,max.skewness,length.out=n)

    tmp <- list(x=b1seq,y=tnofn(b1seq^2))                                       # type I
    if ((max(tmp$y)<max.kurtosis)&&(t3fn(max.skewness^2)>max.kurtosis)) {
      tmp$x <- c(tmp$x,max.skewness^2)
      tmp$y <- c(tmp$y,max.kurtosis)
    }
    tmp$x <- c(tmp$x,rev(b1seq)); tmp$y <- c(tmp$y,t3fn(rev(b1seq^2)))
    polygon(tmp,col="yellow")

    tmp <- list(x = b1seq[b1seq<sqrt(32)], y = t5fn(b1seq[b1seq<sqrt(32)]^2))   # type IV
    if (max(tmp$y)<max.kurtosis) {
      tmp$x <- c(tmp$x,max.skewness^2)
      tmp$y <- c(tmp$y,max.kurtosis)
    }
    tmp$x <- c(tmp$x,0) 
    tmp$y <- c(tmp$y,max.kurtosis)
    polygon(tmp,col="lightcyan")

    tmp <- list(x=b1seq,y=t3fn(b1seq^2))                                        # type VI
    if ((max(tmp$y)<max.kurtosis)&&(t5fn(max.skewness^2)>max.kurtosis)) {
      tmp$x <- c(tmp$x,max.skewness^2)
      tmp$y <- c(tmp$y,max.kurtosis)
    }
    tmp$x <- c(tmp$x, rev(b1seq[b1seq<sqrt(32)]))
    tmp$y <- c(tmp$y, t5fn(rev(b1seq[b1seq<sqrt(32)]^2)))
    polygon(tmp,col="pink")

    curve(tnofn(x^2),from=0,to=max.skewness,n=n,add=TRUE,col="black",lwd=lwd)   # no Pearson
    curve(t3fn(x^2),from=0,to=max.skewness,n=n,add=TRUE,col="red",lwd=lwd)      # type III
    curve(t5fn(x^2),from=0,to=min(sqrt(31.999),max.skewness^2),n=n,add=TRUE,    # type V 
          col="blue", lwd=lwd)  
    lines(c(0,0),c(1,3),col="orange",lwd=lwd,xpd=NA)                            # type II
    lines(c(0,0),c(3,max.kurtosis),col="cyan",lwd=lwd,xpd=NA)                   # type VII
    points(0,3,pch=16,col="brown")                                              # type 0
  }
  if (legend) legend("bottomright",
                     c("2-point-distr.","Pearson 0","Pearson I","Pearson II",
                       "Pearson III","Pearson IV","Pearson V","Pearson VI",
                       "Pearson VII"),
                     col=c("black","brown","yellow","orange","red","lightcyan",
                           "blue","pink","cyan"),
                     lwd=c(lwd,NA,NA,lwd,lwd,NA,lwd,NA,lwd),
                     lty=c(1,NA,NA,1,1,NA,1,NA,1),
                     pch=c(NA,16,15,NA,NA,15,NA,15,NA),
                     pt.cex=c(1,1,2,1,1,2,1,2,1))
}
