norm.OC <- function(k = NULL, alpha = NULL, P = NULL, n, side = 1, 
                    method = c("HE", "HE2", "WBE", "ELL", "KM", "EXACT", "OCT"), m = 50){
  if(side != 1 && side != 2){
    stop(paste("Must specify a one-sided or two-sided procedure!", 
               "\n"))
  }  
  if(length(n)<2){
    stop(paste("'n' needs to be a vector of at least length 2 to produce an OC curve.", 
               "\n"))
  }  
  n <- sort(n)
  method <- match.arg(method)  
  col.blind=c("#000000", "#E69F00", "#56B4E9", "#009E73", 
                  "#F0E442", "#0072B2", "#D55E00", "#CC79A7",
                  "#7FFF00","#7D26CD")
  if(is.null(P)){
    if(length(k)!=1|length(alpha)<1){
      stop(paste("Check values specified for k, n, and alpha!", 
                 "\n"))
    }
    if(length(alpha)>10){
      warning("Too many values of alpha specified!  Using only the first 10 values.",call.=FALSE)
    }
    dev.new(width=11,height=5)
    par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
    alpha <- sort(alpha)[1:min(length(alpha),10)]
    tmp.obj = paste("(k=",k,")",sep="")
    all.P <- cbind(sapply(1:length(alpha), function(i) sapply(1:length(n), function(j) uniroot(function(P) k-K.factor(n=n[j],alpha=alpha[i],P=P,method=method,m=m,side=side),lower=1e-10,upper=1-1e-10)$root)))
    plot(n,rep(0,length(n)),ylab="P",xlab="n",col=0,main=bquote("Normal Tolerance Interval OC Curve for P" ~ .(tmp.obj)),ylim=c(min(all.P),1))
    rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = colors()[356])
    for(i in 1:length(alpha)) lines(n,all.P[,i],type="o",pch=19,col=col.blind[i],cex=.5)
    legend("topright",legend=formatC(1-alpha,width=5,format='f',digits=3,flag='0')[length(alpha):1],title=expression(paste("     1-",alpha)),text.col=col.blind[length(alpha):1],bty="n",title.col="navy",inset=c(-0.125,0))
  } else if(is.null(alpha)){
    if(length(k)!=1|length(P)<1){
      stop(paste("Check values specified for k, n, and P!", 
                 "\n"))
    }
    if(length(P)>10){
      warning("Too many values of P specified!  Using only the first 10 values.",call.=FALSE)
    }    
    dev.new(width=11,height=5)
    par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
    P <- sort(P)[1:min(length(P),10)]
    tmp.obj = paste("(k=",k,")",sep="")
    all.alpha <- sapply(1:length(P), function(i) sapply(1:length(n), function(j) uniroot(function(alpha) k-K.factor(n=n[j],alpha=alpha,P=P[i],method=method,m=m,side=side),lower=1e-10,upper=1-1e-10)$root))
    plot(n,rep(0,length(n)),ylab=expression(paste("1-",alpha,sep="")),xlab="n",col=0,main=bquote("Normal Tolerance Interval OC Curve for 1-" ~ alpha ~ .(tmp.obj)),ylim=c(min(1-all.alpha),1))
    rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = colors()[356])
    for(i in 1:length(P)) lines(n,1-all.alpha[,i],type="o",pch=19,col=col.blind[i],cex=.5)
    legend("topright",legend=formatC(P,width=5,format='f',digits=3,flag='0'),title="     P",text.col=col.blind,bty="n",title.col="navy",inset=c(-0.125,0))
  } else if(is.null(k)){
    if((length(P)*length(alpha))>10){
      warning("Too many combinations of alpha and P specified!  Using only the first 10 such combinations.",call.=FALSE)
    }    
    dev.new(width=11,height=5)
    par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
    alpha <- sort(alpha)
    P <- sort(P)
    tmp <- K.table(n=n,alpha=alpha,P=P,method=method,m=m,side=side)
    all.k <- sapply(tmp,c)
    all.k <- all.k[1:min(nrow(all.k),10),]
    tmp.alpha <- rep(1-alpha,length(P))
    tmp.P <- sort(rep(P,length(alpha)))
    labs <- cbind(formatC(tmp.alpha,width=5,format='f',digits=3,flag='0'),formatC(tmp.P,width=5,format='f',digits=3,flag='0'))[1:min(nrow(all.k),10),]
    new.labs <- NULL
    for(i in 1:nrow(labs)) new.labs <- c(new.labs,paste("(",labs[i,1],",",labs[i,2],")",sep=""))
    plot(n,rep(0,length(n)),ylab="k",xlab="n",col=0,main=bquote("Normal Tolerance Interval OC Curve for k and n" ~ ""),ylim=c(0,max(all.k)))
    rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = colors()[356])
    for(i in 1:nrow(all.k)) lines(n,all.k[i,],type="o",pch=19,col=col.blind[i],cex=.5)
    legend("topright",legend=new.labs,title=expression(paste("     (1-",alpha,",P)",sep="")),text.col=col.blind,bty="n",title.col="navy",inset=c(-0.15,0))
  } else{
    stop(paste("Check values specified for k, n, alpha, and P!", 
               "\n"))
  }
}
