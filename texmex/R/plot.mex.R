plot.mex <- function(x,quantiles=seq(0.1,by=0.2,len=5),col="grey",...){
   if ( class( x ) != "mex" )
         stop( "you need to use an object with class 'mex'" )

   mar <- x[[1]]
   dep <- x[[2]]
   z <- dep$Z
   n <- dim(z)[1]
   xmax <- max(mar$data[,dep$which])
   sig <- coef(mar)[3,dep$which]
   xi <- coef(mar)[4,dep$which]
   marThr <- mar$mth[dep$which]
   marP   <- mar$mqu[dep$which]
   if(xi < 0) upperEnd <- marThr - sig/xi
   len <- 1001
   
   for(i in 1:(dim(z)[2])){
      p <- seq(dep$dqu,1-1/n,length=n)
      plot(p,z[,i],xlab=paste("F(",dep$conditioningVariable,")",sep=""),
           ylab=paste("Z   ",colnames(z)[i]," | ",dep$conditioningVariable,sep=""),col=col,...)
      lines(lowess(p,z[,i]),col=2)
      plot(p,abs(z[,i] - mean(z[,i])),xlab=paste("F(",dep$conditioningVariable,")",sep=""),
           ylab=paste("|Z - mean(Z)|   ",colnames(z)[i]," | ",dep$conditioningVariable,sep=""),col=col,...)
      lines(lowess(p,abs(z[,i] - mean(z[,i]))),col=4)

      depThr <- c(quantile(mar$data[,dep$which],dep$dqu))
      d <- xmax-depThr
      xlim <- c(depThr-0.1*d, depThr + 1.5*d)
      names(xlim) <- NULL

      SetPlim <- TRUE
      if(xi < 0 && xlim[2] > upperEnd){
        xlim[2] <-  upperEnd
        plim <- 1
        SetPlim <- FALSE
      }
      
      if( SetPlim ) plim <- pgpd(xlim[2],sigma=sig,xi=xi,u=marThr)
      
      plotp <- seq(dep$dqu,plim,len=len)[-len] # take out largest point to avoid Inf in next line if xlim[2]==upperEnd
      co <- coef(dep)[,i]
      xq <- -log(-log(plotp))
      zq <- quantile(dep$Z[,i],quantiles)
      yq <- sapply(zq, function(z, co, xq){
      						co["a"] * xq + co["c"] - co["d"]*log(xq) + xq^co["b"] * z
      					}, # Close function
      					 xq=xq, co=co
      			   ) # Close sapply
      
      plotx <- revTransform(plotp,data=mar$data[,dep$which],qu=marP,th=marThr,sigma=sig,xi=xi)
      ploty <- apply(exp(-exp(-yq)),2,revTransform,data=as.matrix(mar$data[,-dep$which])[,i],
                     qu=mar$mqu[-dep$which][i],th=mar$mth[-dep$which][i],
                     sigma=coef(mar)[3,-dep$which][i],xi=coef(mar)[4,-dep$which][i])
      
      plot(mar$data[,dep$which],as.matrix(mar$data[,-dep$which])[,i], xlab=dep$conditioningVariable,ylab=colnames(z)[i],col=col,...)
      abline(v=depThr)
      for(j in 1:length(quantiles)){
        lines(plotx,ploty[,j],lty=2)
      }
    }
}

  
