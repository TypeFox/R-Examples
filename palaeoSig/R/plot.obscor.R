#plot.obscor <-
#function(x, xlab="WA optima", ylab="RDA scores", f=1,...){
#  plot(x=x$x$Optima,y=x$x$RDA1, cex=x$x$abun*f, xlab=xlab, ylab=ylab,...)
#}



plot.obscor<-function (x, xlab, ylab, f=5, which=1, label="env",abun="abun.calib",p.val=0.95,...) 
{
    weightings <- c("abun.fos", "abun.calib", "abun.joint", "n2.fos", "n2.calib", "n2.joint", "unweighted")
    w <- pmatch(abun, weightings)
    if (is.na(w)) {
        stop("Unknown abundance weighting")
    }
    w<-weightings[w]


  if(which==1){
    if(missing(xlab))xlab= "WA optima"
    if(missing(ylab))ylab = "RDA scores"            
    if(w=="unweighted")a<-rep(1, nrow(x$ob$x))
    else{
      a<-x$ob$x[[w]]
      a<-a/max(a)*f
    }
    plot(x = x$ob$x$Optima, y = x$ob$x$RDA1, cex = a, xlab = xlab, ylab = ylab, ...)
  }else if(which==2){
    if(missing(xlab))xlab = ifelse(w=="unweighted","Correlation","Weighted correlation")
    sim<-x$sim[[w]]
    ob<-x$ob$res[w]
    hist(sim, xlim = range(c(sim,ob)), xlab = xlab, col = "grey80", border = NA, ...)
    abline(v = ob, col = 1)
    abline(v = quantile(sim, prob=p.val), col = 2, lty=3)
    
    text(ob, par()$usr[4] * 0.9, label = label, srt = 90, pos = 2)
    box()
  }else stop("which==what")
}
