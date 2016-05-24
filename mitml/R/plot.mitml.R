plot.mitml <- function(x, print=c("beta","psi","sigma"), group="all", trace=c("imputation","burnin","all"), smooth=3, export=c("none","png","pdf"), dev.args=list(), ...){
# plot method for objects of class "mitml"

  # .movingAverage <- mitml:::.movingAverage

  # retrieve data and variable names
  vrs <- x$model
  clus <- x$model$clus
  pvrs <- seq_along(attr(vrs,"full.names")$pvrs)
  qvrs <- seq_along(attr(vrs,"full.names")$qvrs)
  names(pvrs) <- attr(vrs,"full.names")$pvrs
  names(qvrs) <- attr(vrs,"full.names")$qvrs

  # check for random L1
  rl1 <- x$random.L1=="full"

  prt <- match.arg(print,several.ok=TRUE)
  shw <- match.arg(trace)
  export <- match.arg(export)

  # grouping
  grp.labels <- unique(attr(x$data,"group"))
  grp.imp <- length(grp.labels)
  if(class(group)=="numeric") grp.labels <- grp.labels[group]
  grp <- length(grp.labels)

  # export, graphical parameters
  if(export!="none"){
    wd <- getwd()
    out <- file.path(wd,"mitmlPlots")
    if(!file.exists(out)) dir.create(out)
  }else{
    do.call(dev.new,dev.args)
    devAskNewPage(ask=FALSE)
  }
  oldpar <- par(no.readonly=TRUE)

  # ***
  # start plotting
  #

  for(gg in 1:grp){

  # grouping
  if(grp.imp>1){
    glab <- paste(",Group:",grp.labels[gg],sep="")
    gfile <- paste("Group-",grp.labels[gg],"_",sep="")
  }else{
    glab <- gfile <- ""
  }

  # expand yvrs for multiple categories
  yvrs <- vrs$yvrs
  cvrs <- attr(x$data,"cvrs")
  # expand categorical target variables
  if(!is.null(cvrs)){
    yvrs <- c(yvrs[!yvrs%in%cvrs], cvrs)
    for(cc in 1:length(cvrs)){
      cv <- cvrs[cc]
      ci <- which(yvrs==cv)
      yi <- 1:length(yvrs)
      nlev <- attr(x$data,"levels")[gg,cc]
      if(nlev>2){
        newy <- paste0(cv, 1:(nlev-1))
      }else{
        newy <- cv
      }
      sel0 <- yi[yi<ci]
      sel1 <- yi[yi>ci]
      yvrs <- c(yvrs[sel0],newy,yvrs[sel1])
    }
  }
  ynam <- yvrs
  yvrs <- seq_along(yvrs)
  names(yvrs) <- ynam

  # plots for fixed regression coefficients
  if("beta" %in% prt){
  for(ii in yvrs){
    for(jj in pvrs){
    
      if(export!="none"){
        filename <- paste("BETA_",gfile,names(yvrs[ii]),"_ON_",names(pvrs[jj]),".",export,sep="")
        filename <- gsub("[(),]","",filename)
        filename <- gsub("[[:space:]]","-",filename)
        out.args <- c(list(file=file.path(out,filename)),dev.args)
        do.call(export, out.args)
      }
      
      layout(matrix(c(1,2,3,4),2,2), c(5,1), c(1.13,1))

      switch(shw,
        imputation={trc <- x$par.imputation[["beta"]][jj,ii,,gg]},
        burnin={trc <- x$par.burnin[["beta"]][jj,ii,,gg]},
        all={trc <- c(x$par.burnin[["beta"]][jj,ii,,gg], x$par.imputation[["beta"]][jj,ii,,gg])}
      )

      # trace plot
      par(mar=c(3,3,2,0)+0.5, mgp=c(2,1,0), font.lab=2)
      ymin <- min(trc)
      ymax <- max(trc)
      yr <- ymax-ymin
      plot(trc, type="l", ylab="Trace", xlab="Iteration", xaxt="n", ylim=c(ymin-yr*.03, ymax+yr*.03),
           ...)
      axt <- axTicks(1)
      title(main=paste("Beta [",jj,",",ii,glab,"]: ",names(yvrs[ii])," ON ",names(pvrs[jj]),sep=""), cex.main=1)
      if(shw=="imputation"){
        axl <- sprintf("%d", axt + dim(x$par.burnin[["beta"]])[3])
      }else{
        axl <- sprintf("%d", axt)
      }
      axis(side=1, at=axt, labels=axl)
      if(shw=="all") abline(v=dim(x$par.burnin[["beta"]])[3], col="blue")

      # trend line for trace (moving window average)
      if(all(is.numeric(smooth),smooth>0)){
        B <- floor(x$iter$iter/smooth)
        mwa <- .movingAverage(trc,B,fill=TRUE)
        lines(mwa, col="grey65")
      }

      # autocorrelation plot
      par(mar=c(3,3,1,0)+0.5)
      drw <- x$par.imputation[["beta"]][jj,ii,,gg]
      acf(drw, lag.max=x$iter[["iter"]], ylim=c(-.1,1), yaxt="n", main=NULL, ylab="ACF", ci=0, ...)
      axis(side=2, at=c(0,.5,1))
      abline(h=c(-.1,.1), col="blue")

      # kernel density plot
      par(mar=c(3,0,2,0)+0.5, mgp=c(2,0,0))
      ddrw <- density(drw)
      plot(x=ddrw$y, y=ddrw$x, type="l", xaxt="n", yaxt="n", xlab="", ylab="", ylim=c(ymin-yr*.03, ymax+yr*.03), ...)

      # posterior summary
      par(mar=c(3,-0.5,2,-0.5)+0.5)
      plot.new()
      text(0,0.5,paste("EAP:   ", sprintf(fmt="%.3f", mean(drw)), "\n",
                       "MAP:   ", sprintf(fmt="%.3f", ddrw$x[which.max(ddrw$y)]), "\n",
                       "SD:    ", sprintf(fmt="%.3f", sd(drw)), "\n",
                       "2.5%:  ", sprintf(fmt="%.3f", quantile(drw,.025)), "\n",
                       "50%:   ", sprintf(fmt="%.3f", quantile(drw,.5)), "\n",
                       "97.5%: ", sprintf(fmt="%.3f", quantile(drw,.975)), "\n",
                       sep=""), adj=c(0,.5), cex=.8, family="mono", font=2, ...)
    
      if(export!="none"){
        dev.off()
      }else{
        devAskNewPage(ask=TRUE)
      }
  }}}

  # plots for random effects
  if("psi" %in% prt){
  bvec <- t(expand.grid(qvrs, yvrs))
  for(ii in 1:(length(yvrs)*length(qvrs))){
    for(jj in ii:(length(yvrs)*length(qvrs))){

      if(export!="none"){
        filename <- paste("PSI_",gfile,names(yvrs[bvec[2,ii]]),"_ON_",names(qvrs[bvec[1,ii]]),
          "_WITH_",names(yvrs[bvec[2,jj]]),"_ON_",names(qvrs[bvec[1,jj]]),".",export,sep="")
        filename <- gsub("[(),]","",filename)
        filename <- gsub("[[:space:]]","-",filename)
        out.args <- c(list(file=file.path(out,filename)),dev.args)
        do.call(export, out.args)
      }

      layout(matrix(c(1,2,3,4),2,2), c(5,1), c(1.13,1))

      switch(shw,
        imputation={trc <- x$par.imputation[["psi"]][jj,ii,,gg]},
        burnin={trc <- x$par.burnin[["psi"]][jj,ii,,gg]},
        all={trc <- c(x$par.burnin[["psi"]][jj,ii,,gg], x$par.imputation[["psi"]][jj,ii,,gg])}
      )

      # trace plot
      par(mar=c(3,3,2,0)+0.5, mgp=c(2,1,0), font.lab=2)
      ymin <- min(trc)
      ymax <- max(trc)
      yr <- ymax-ymin
      plot(trc, type="l", ylab="Trace", xlab="Iteration", xaxt="n", ylim=c(ymin-yr*.03, ymax+yr*.03), 
           ...)
      title(main=paste("Psi [",ii,",",jj,glab,"]: (",names(yvrs[bvec[2,ii]])," ON ",names(qvrs[bvec[1,ii]]),
            ") WITH (",names(yvrs[bvec[2,jj]])," ON ",names(qvrs[bvec[1,jj]]),")",sep=""), cex.main=1)
      axt <- axTicks(1)
      if(shw=="imputation"){
        axl <- sprintf("%d", axt + dim(x$par.burnin[["psi"]])[3])
      }else{
        axl <- sprintf("%d", axt)
      }
      axis(side=1, at=axt, labels=axl)
      if(shw=="all") abline(v=dim(x$par.burnin[["psi"]])[3], col="blue")

      # trend line for trace (moving window average)
      if(all(is.numeric(smooth),smooth>0)){
        B <- floor(x$iter$iter/smooth)
        mwa <- .movingAverage(trc,B,fill=TRUE)
        lines(mwa, col="grey65")
      }

      # autocorrelation plot
      par(mar=c(3,3,1,0)+0.5)
      drw <- x$par.imputation[["psi"]][jj,ii,,gg]
      acf(drw, lag.max=x$iter[["iter"]], ylim=c(-.1,1), yaxt="n", main=NULL, ylab="ACF", ci=0, ...)
      axis(side=2, at=c(0,.5,1))
      abline(h=c(-.1,.1), col="blue")

      # kernel density plot
      par(mar=c(3,0,2,0)+0.5, mgp=c(2,0,0))
      ddrw <- density(drw)
      plot(x=ddrw$y, y=ddrw$x, type="l", xaxt="n", yaxt="n", xlab="", ylab="", ylim=c(ymin-yr*.03, ymax+yr*.03), ...)

      # posterior summary
      par(mar=c(3,-0.5,2,-0.5)+0.5)
      plot.new()
      text(0,0.5,paste("EAP:   ", sprintf(fmt="%.3f", mean(drw)), "\n",
                       "MAP:   ", sprintf(fmt="%.3f", ddrw$x[which.max(ddrw$y)]), "\n",
                       "SD:    ", sprintf(fmt="%.3f", sd(drw)), "\n",
                       "2.5%:  ", sprintf(fmt="%.3f", quantile(drw,.025)), "\n",
                       "50%:   ", sprintf(fmt="%.3f", quantile(drw,.5)), "\n",
                       "97.5%: ", sprintf(fmt="%.3f", quantile(drw,.975)), "\n",
                       sep=""), adj=c(0,.5), cex=.8, family="mono", font=2, ...)

    if(export!="none"){
      dev.off()
    }else{
      devAskNewPage(ask=TRUE)
    }
  }}}

  # plots for residuals
  if("sigma" %in% prt){

  clus2 <- unique(x$data[,clus])
  clus3 <- if(rl1) seq_along(clus2) else 1
  for(icl in clus3){

  for(ii in 1:length(yvrs)){
    for(jj in ii:length(yvrs)){

      if(export!="none"){
        filename <- paste0("SIGMA_",gfile,names(yvrs[ii]),"_WITH_",names(yvrs[jj]), 
                           if(rl1) paste0(clus,clus2[icl]),".",export)
        filename <- gsub("[(),]","",filename)
        filename <- gsub("[[:space:]]","-",filename)
        out.args <- c(list(file=file.path(out,filename)),dev.args)
        do.call(export, out.args)
      }

      layout(matrix(c(1,2,3,4),2,2), c(5,1), c(1.13,1))

      jj2 <- jj+(icl-1)*length(yvrs)
      switch(shw,
        imputation={
          trc <- x$par.imputation[["sigma"]][jj2,ii,,gg]
        },
        burnin={
          trc <- x$par.burnin[["sigma"]][jj2,ii,,gg]
        },
        all={
          trc <- c(x$par.burnin[["sigma"]][jj2,ii,,gg], x$par.imputation[["sigma"]][jj2,ii,,gg])
        }
      )

      # trace plots
      par(mar=c(3,3,2,0)+0.5, mgp=c(2,1,0), font.lab=2)
      ymin <- min(trc)
      ymax <- max(trc)
      yr <- ymax-ymin
      plot(trc, type="l", ylab="Trace", xlab="Iteration", xaxt="n", ylim=c(ymin-yr*.03, ymax+yr*.03), 
           ...)
      title(main=paste0("Sigma [",ii,",",jj,glab,"]: ",names(yvrs[ii])," WITH ",names(yvrs[jj]),
                        if(rl1) paste0(" [",clus,":",clus2[icl],"]")), cex.main=1)
      axt <- axTicks(1)
      if(shw=="imputation"){
        axl <- sprintf("%d", axt + dim(x$par.burnin[["sigma"]])[3])
      }else{
        axl <- sprintf("%d", axt)
      }
      axis(side=1, at=axt, labels=axl)
      if(shw=="all") abline(v=dim(x$par.burnin[["sigma"]])[3], col="blue")

      # trend line for trace (moving window average)
      if(all(is.numeric(smooth),smooth>0)){
        B <- floor(x$iter$iter/smooth)
        mwa <- .movingAverage(trc,B,fill=TRUE)
        lines(mwa, col="grey65")
      }

      # autocorrelation plot
      par(mar=c(3,3,1,0)+0.5)
      drw <- x$par.imputation[["sigma"]][jj2,ii,,gg]
      acf(drw, lag.max=x$iter[["iter"]], ylim=c(-.1,1), yaxt="n", main=NULL, ylab="ACF", ci=0, ...)
      axis(side=2, at=c(0,.5,1))
      abline(h=c(-.1,.1), col="blue")

      # kernel density plot
      par(mar=c(3,0,2,0)+0.5, mgp=c(2,0,0))
      ddrw <- density(drw)
      plot(x=ddrw$y, y=ddrw$x, type="l", xaxt="n", yaxt="n", xlab="", ylab="", ylim=c(ymin-yr*.03, ymax+yr*.03), ...)

      # posterior summary
      par(mar=c(3,-0.5,2,-0.5)+0.5)
      plot.new()
      text(0,0.5,paste("EAP:   ", sprintf(fmt="%.3f", mean(drw)), "\n",
                       "MAP:   ", sprintf(fmt="%.3f", ddrw$x[which.max(ddrw$y)]), "\n",
                       "SD:    ", sprintf(fmt="%.3f", sd(drw)), "\n",
                       "2.5%:  ", sprintf(fmt="%.3f", quantile(drw,.025)), "\n",
                       "50%:   ", sprintf(fmt="%.3f", quantile(drw,.5)), "\n",
                       "97.5%: ", sprintf(fmt="%.3f", quantile(drw,.975)), "\n",
                       sep=""), adj=c(0,.5), cex=.8, family="mono", font=2, ...)

    if(export!="none"){
      dev.off()
    }else{
      devAskNewPage(ask=TRUE)
    }
  }}}}

  }

  plot.new()
  par(oldpar)
  if(export=="none") devAskNewPage(ask=FALSE)
  dev.off()
  invisible()
}

# moving window average for time series
.movingAverage <- function(x, B, fill=TRUE){

    x1 <- cumsum(x)
    N <- length(x)
    y <- rep(NA,N)
    i <- seq(B+1 , N-B)
    xdiff <- x1[ -seq(1,B) ] - x1[ -seq(N-B+1,N) ] 
    xdiff <- xdiff[ - seq(1,B) ] 

    y[i]  <- ( x1[i] + xdiff - c(0,x1[ -seq(N-2*B,N) ]) ) / (2*B+1)

  # fill NAs at beginning and end of time series
  if(fill){
    j <- seq(0,B-1)
    ybeg <- sapply(j, function(z) sum( x[ seq(1,(2*z+1)) ]) / (2*z+1) )
    yend <- sapply(rev(j), function(z) sum( x[ seq(N-2*z,N) ] ) / (2*z+1) )
    y[j+1] <- ybeg
    y[rev(N-j)] <- yend
  }

  y
}

