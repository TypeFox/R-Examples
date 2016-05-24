centralPlot  <- function(x, p=0.5,col.c='red',col.e='slategray',lty=c(1,3),gradient=FALSE,gradient.ramp=NULL,main=NULL,cex=1,...)
{
    par(mar=c(4,5,3,3),xpd=FALSE)
    p=1-p
    x<-as.matrix(x)
    n <- nrow(x)
    I <- MBD(x, plotting = FALSE)$ordering
    N <- n - floor(p * n)  ## number of central curves
    if (ncol(x) == 1)  x <- t(x)
    if ((p<0)|(p>1)) stop("Incorrect proportion of curves to discard...")
    if (length(col.c)<N) col.c<-rep(col.c,length.out=N)
    if (length(col.e)<(n-N)) col.e<-rep(col.e,length.out=n-N)
    if (gradient) {
         if (is.null(gradient.ramp)) {gradient.ramp <- colorRampPalette(c('gold','red'))(N)}  
         else {
           if (length(gradient.ramp)<2) gradient.ramp<-rep(gradient.ramp,N)  ## one color, N times
           else gradient.ramp<-colorRampPalette(c(gradient.ramp[2],gradient.ramp[1]))(N)
           }  ## end ELSE
         col.c <- gradient.ramp
         }  ## end IF gradient
    if (length(lty)<2) stop ("Not enough line types defined...")
    if (is.null(main)) main<-'Central plot'
    if (N>0) {
      m1 <- as.matrix(x[I[N:1], ])
      if (ncol(m1)==1) m1 <-t(m1)
      }  ## end IF
    else {
      m1<-matrix()
      warning("Too small p: all samples discarded...")
      } ## end ELSE
    if (N<n){
      m2 <- x[I[(N+1):n], ]
      Gene.Expression<-t(m2)
      matplot(Gene.Expression,type="l", lty=lty[2],col=col.e,xlab='',ylim=c(min(x),max(x)),main=main,...)
      if (N>0) {
         matlines(t(m1),lty=lty[1],col=col.c,...)
         }
      } 
    else {
      Gene.Expression<-t(m1)
      matplot(Gene.Expression,type="l", lty=lty[1],col=rev(col.c),xlab='',ylim=c(min(x),max(x)),main=main)
      warning("No external curves discarded...")
      }  ## end ELSE
    par(xpd=TRUE)
    legend("top",lty=lty,col=c(col.c[N],col.e),legend=c('deepest curve(s)',paste(p*100,'% most external curves')),cex=cex)

    }
