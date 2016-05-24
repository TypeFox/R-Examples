plot.hclustvar <-
  function(x,type="tree",which=c(1:2), ask = prod(par("mfcol")) < length(which) && dev.interactive(),sub="", ...){
    if (!inherits(x, "hclustvar")) 
      stop("use only with \"hclustvar\" objects")
    if (!is.numeric(which) || any(which < 1) || any(which > 2)) 
      stop("'which' must be in 1:2")
    class(x) <- "hclust"
    #show <- rep(FALSE, 2)
    #show[which] <- TRUE
    #if (ask) {
    #  oask <- devAskNewPage(TRUE)
    #  on.exit(devAskNewPage(oask))
    # }	
    #if (show[1]) plot(x, hang=-1, xlab="", sub=sub, ...)	
    #if (show[2]) {
    #	plot(x=seq(length(x$height),1),x$height,xaxt = "n",ylim=c(0,max(x$height)),xlab="number of clusters",ylab="Height",main="Aggregation levels",type="n")
    #	points(x=seq(length(x$height),1),x$height,pch=3)
    #	axis(side=1,at=seq(1,length(x$height)),labels=paste(1:(length(x$height))))
    #}
    if ((type !="tree") && (type !="index")) { stop("type must be equal to  \"tree \" or  \"index\"")} else {
      if (type=="tree") plot(x, hang=-1, xlab="", sub=sub, ...) 
      if (type=="index"){
        plot(x=seq(length(x$height),1),x$height,xaxt = "n",ylim=c(0,max(x$height)),xlab="number of clusters",ylab="Height",main="Aggregation levels",type="n")
        points(x=seq(length(x$height),1),x$height,pch=3)
        axis(side=1,at=seq(1,length(x$height)),labels=paste(1:(length(x$height)))) 
       }
    }
    
    class(x) <- c("hclustvar","hclust")
  }
