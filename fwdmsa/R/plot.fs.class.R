plot.fs.class <-
function(x, 
         type = "objective", 
         observations = all.observations, 
         id.observation = FALSE, 
         items = all.items, 
         id.item = FALSE,
         step = default.step,
         reference.step=default.reference.step,
         id.scale = default.scale,
         tukey.fences = TRUE,
         add = FALSE,
         n0 = FALSE,
         n1 = FALSE,
         n2 = FALSE,
         lower.c = default.lower.c,
         col = default.col,
         lwd = default.lwd,
         lty=default.lty,
         ylim=default.ylim,
         xlim=default.xlim, ...){
# x is output from fs.MSA function

  # READ 
   N <- x$data$N
   J <- x$data$J
   m <- x$data$m

  # DETERMINE observations and MINSIZE
   all.observations <- 1:N
   all.items <- 1:J
  
 # DETERMINE PLOT
#   if (substr(type,1,1)=="O" || substr(type,1,1)=="o") type <- "objective" else
   if (substr(type,1,2)=="Mi"|| substr(type,1,2)=="mi") type <- "minexcl"   else
   if (substr(type,1,2)=="Ma"|| substr(type,1,2)=="ma") type <- "maxincl"   else
   if (substr(type,1,1)=="G"|| substr(type,1,1)=="g")   type <- "gap"       else
   if (substr(type,1,1)=="H"|| substr(type,1,1)=="h" ||substr(type,1,1)=="C"|| substr(type,1,1)=="c")   type <- "coefH"     else
   if (substr(type,1,1)=="R"|| substr(type,1,1)=="r")   type <- "restscore" else
   if (substr(type,1,1)=="I"|| substr(type,1,1)=="i")   type <- "IRF"       else 
   if (substr(type,1,1)=="F"|| substr(type,1,1)=="f")   type <- "followup"  else
   if (substr(type,1,1)=="N"|| substr(type,1,1)=="n")   type <- "num.scale"  else
   if (substr(type,1,1)=="S"|| substr(type,1,1)=="s")   type <- "scale"  else type <- "objective"
# PLOT "objective" 
   if (type=="objective"){
     observation.obj <- x$objective$observation.obj
       
     if(add==FALSE){
      if(is.numeric(id.observation)==FALSE){
       default.ylim <- c(0,max(observation.obj,na.rm=TRUE))
       default.xlim <- c(0,N)
       default.col <- 1
       default.lty <- 1
       default.lwd <- 1
       plot(NA, xlim=xlim, ylim=ylim, ylab="Objective function", xlab="Subsample size n",las=1)
       for(i in observations){lines(observation.obj[i,], col=col, lty=lty, lwd=lwd)}
      }     
     
      if(is.numeric(id.observation)==TRUE){
       default.ylim <- c(0,max(observation.obj,na.rm=TRUE))
       default.xlim <- c(0,N+N/20)
       default.col <- 2
       default.lty <- 1
       default.lwd <- 3
       plot(NA, xlim=xlim, ylim=ylim, ylab="Objective function", xlab="Subsample size n",las=1)
       for(i in observations){lines(observation.obj[i,], lty=i)}
       for(i in 1:length(id.observation)){
        lines(observation.obj[id.observation[i],], lwd=lwd, col=col, lty=lty)
        text(N,observation.obj[id.observation[i],N], id.observation[i], cex=1, pos=4)
       }
      }
    if(n0==TRUE){axis(1, at=x$initial.subsample.size, labels=expression(italic(n[0])), tick=TRUE, padj=-1.3, cex.axis=1)}
    if(is.numeric(n1)==TRUE){axis(1, at=n1, labels=expression(italic(n[1])), tick=TRUE, padj=-1.3, cex.axis=1)}
    if(n2==TRUE){axis(1, at=x$n2, labels=expression(italic(n[2])), tick=TRUE, padj=-1.3, cex.axis=1)}
    }
    
    if(add==TRUE){
      if(is.numeric(id.observation)==FALSE){
       default.ylim <- c(0,max(observation.obj,na.rm=TRUE))
       default.xlim <- c(0,N)
       default.col <- 1
       default.lty <- 1
       default.lwd <- 1
       for(i in observations){lines(observation.obj[i,], col=col, lty=lty, lwd=lwd)}
      }     
     
      if(is.numeric(id.observation)==TRUE){
       default.ylim <- c(0,max(observation.obj,na.rm=TRUE))
       default.xlim <- c(0,N+N/20)
       default.col <- 2
       default.lty <- 1
       default.lwd <- 3
       for(i in 1:length(id.observation)){
        lines(observation.obj[id.observation[i],], lwd=lwd, col=col, lty=lty)
        text(N,observation.obj[id.observation[i],N], id.observation[i], cex=1, pos=4)
       }
      }
    }
   }


# PLOT "minexcl" 
   if (type=="minexcl"){
    default.ylim <- c(0, max(x$objective$observation.obj,na.rm=TRUE))
    default.xlim <- c(0, N)
    
    min.excl <- x$objective$min.excl
    default.col <- 1
    default.lty <- 1
    default.lwd <- 2
        
    if(tukey.fences==TRUE){
     if(add==FALSE){
      plot(min.excl, ylim=ylim, xlim=xlim, type="l", xlab="Subsample size n", ylab="Objective function", las=1, col=col, lwd=lwd, lty=lty)
      lines(x$tukey.upper.fence[,1],col=2,lwd=2)
      #lines(x$tukey.upper.fence[,2],col=2,lwd=2,lty=2)
    
     }
     if(add==TRUE){
      lines(min.excl, ylim=ylim, xlim=xlim, type="l", xlab="Subsample size n", ylab="Objective function", las=1, col=col, lwd=lwd, lty=lty)
      lines(x$tukey.upper.fence[,1],col=2,lwd=2)
      #lines(x$tukey.upper.fence[,2],col=2,lwd=2,lty=2)
     }
    }
    if(tukey.fences==FALSE){
     if(add==FALSE){
      plot(min.excl, ylim=ylim, xlim=xlim, type="l", xlab="Subsample size n", ylab="Objective function", las=1, col=col, lwd=lwd, lty=lty)
      axis(1, at=x$n2, labels=expression(italic(n[2])), tick=TRUE,padj=-1.3, cex.axis=1)
     }
     if(add==TRUE){
      lines(min.excl, ylim=ylim, xlim=xlim, type="l", xlab="Subsample size n", ylab="Objective function", las=1, col=col, lwd=lwd, lty=lty)
     }
    }
    if(n0==TRUE){axis(1, at=x$initial.subsample.size, labels=expression(italic(n[0])), tick=TRUE, padj=-1.3, cex.axis=1)}
    if(is.numeric(n1)==TRUE){axis(1, at=n1, labels=expression(italic(n[1])), tick=TRUE, padj=-1.3, cex.axis=1)}
    if(n2==TRUE){axis(1, at=x$n2, labels=expression(italic(n[2])), tick=TRUE, padj=-1.3, cex.axis=1)}
   }
   
# PLOT "maxincl" 
   if (type=="maxincl"){
    default.ylim <- c(0, max(x$objective$observation.obj,na.rm=TRUE))
    default.xlim <- c(0, N)
    
    max.incl <- x$objective$max.incl
    default.col <- 1
    default.lty <- 1
    default.lwd <- 2
       
    if(add==FALSE){
     plot(max.incl, ylim=ylim, xlim=xlim, type="l", xlab="Subsample size n", ylab="Objective function", las=1,col=col,lwd=lwd,lty=lty)
     if(n0==TRUE){axis(1, at=x$initial.subsample.size, labels=expression(italic(n[0])), tick=TRUE, padj=-1.3, cex.axis=1)}
     if(is.numeric(n1)==TRUE){axis(1, at=n1, labels=expression(italic(n[1])), tick=TRUE, padj=-1.3, cex.axis=1)}
     if(n2==TRUE){axis(1, at=x$n2, labels=expression(italic(n[2])), tick=TRUE, padj=-1.3, cex.axis=1)}
    }
    
    if(add==TRUE){
     lines(max.incl, ylim=ylim, xlim=xlim, type="l", xlab="Subsample size n", ylab="Objective function", las=1,col=col,lwd=lwd,lty=lty)
    }
   }

# PLOT "gap" 
   if (type=="gap"){
    default.ylim <- c(min(0,min(x$objective$min.excl-x$objective$max.incl,na.rm=TRUE)),max(x$objective$min.excl-x$objective$max.incl,na.rm=TRUE)) 
    default.xlim <- c(0, N)
    default.col <- 1
    default.lty <- 1
    default.lwd <- 2
    
    gap <- x$objective$min.excl-x$objective$max.incl
       
     if(add==FALSE){
       plot(gap, ylim=ylim, xlim=xlim, type="l",las=1, xlab="Subsample size n", ylab="Gap", lwd=lwd, col=col, lty=lty)
       if(n0==TRUE){axis(1, at=x$initial.subsample.size, labels=expression(italic(n[0])), tick=TRUE, padj=-1.3, cex.axis=1)}
       if(is.numeric(n1)==TRUE){axis(1, at=n1, labels=expression(italic(n[1])), tick=TRUE, padj=-1.3, cex.axis=1)}
       if(n2==TRUE){axis(1, at=x$n2, labels=expression(italic(n[2])), tick=TRUE, padj=-1.3, cex.axis=1)}
       }
     if(add==TRUE){
       lines(gap, ylim=ylim, xlim=xlim, type="l", las=1, xlab="Subsample size n", ylab="Gap", lwd=lwd, col=col, lty=lty)
       }
     segments(0, 0, N, 0, col=1, lwd=1)
   }

# PLOT "coefH" 
    if (type=="coefH"){

    default.lower.c <- .3
    default.ylim <- c(0,1)
    default.col <- 1
    default.lty <- 1
    default.lwd <- 2
    
    summary.monotonicity.violation=matrix(,J,N)
    for(i in x$initial.subsample.size:N){summary.monotonicity.violation[,i] <- x$monotonicity[[i]]$Hi}
    if(add==FALSE){
    if(id.item==FALSE){
     default.xlim <- c(0,N)
     plot(NA, ylim=ylim, xlim=xlim, xlab="Subsample size n", ylab="", las=1)
     for(j in items){lines(summary.monotonicity.violation[j,], col=j, lty=j, lwd=lwd)}
    }
        
    if(id.item==TRUE){
     default.xlim <- c(0, N+N/15)
     plot(NA, ylim=ylim, xlim=xlim, xlab="Subsample size n", ylab="", las=1)
     for(j in items){lines(summary.monotonicity.violation[j,], col=j, lty=j, lwd=lwd)}
     for(j in items){text(N, summary.monotonicity.violation[j,N], x$monotonicity[[N]]$I.labels[j], pos=4)}
    }
    
    lines(x$statistics[2,], col=col, lty=lty, lwd=2*lwd)
    segments(0, lower.c, N, lower.c, col=1, lwd=lwd)
    title(ylab=expression(italic(H[j])), line=2.5)
    if(n0==TRUE){axis(1, at=x$initial.subsample.size, labels=expression(italic(n[0])), tick=TRUE, padj=-1.3, cex.axis=1)}
    if(is.numeric(n1)==TRUE){axis(1, at=n1, labels=expression(italic(n[1])), tick=TRUE, padj=-1.3, cex.axis=1)}
    if(n2==TRUE){axis(1, at=x$n2, labels=expression(italic(n[2])), tick=TRUE, padj=-1.3, cex.axis=1)}
    }
    
    if(add==TRUE){
    for(j in items){lines(summary.monotonicity.violation[j,], col=j, lty=j, lwd=lwd)}
    if(id.item==TRUE) {for(j in items){text(N,summary.monotonicity.violation[j,N], x$monotonicity[[N]]$I.labels[j], pos=4)}}
    lines(x$statistics[2,], col=col, lty=lty, lwd=2*lwd)
    }
    }


# PLOT "restscore" 
    if (type=="restscore"){
    
    default.ylim <- c(0,m)
    default.col <- 1
    default.lty <- 1
    default.lwd <- 2
    default.step <- N   
    
    if(add==FALSE){
     for(j in items){
      devAskNewPage(ask = TRUE)
      plot(x$monotonicity[[step]]$results[[j]][[2]][,10], type="b", ylim=ylim, xlim=c(1,dim(x$monotonicity[[N]]$results[[j]][[2]])[1]), pch=19, ylab="", xlab="", las=1, lty=lty, lwd=lwd, col=col, cex=1, xaxt="n", main=x$monotonicity[[N]]$I.labels[j])
      title(ylab=expression(italic(paste(bar(X)[j]," | ",R[(j)],sep=" "))), cex.lab=1.2, line=2)                                                            # NEW: R[(j)] ipv R[(-j)] #
      title(xlab=paste("Rest score group item",x$monotonicity[[N]]$results[[j]][1]), line=3)
      for(l in 1:dim(x$monotonicity[[N]]$results[[j]][[2]])[1]){
       axis(1, at=l, labels=paste(x$monotonicity[[N]]$results[[j]][[2]][l,2], "-", x$monotonicity[[N]]$results[[j]][[2]][l,3], sep=""), tick=TRUE, padj=-.5)
       axis(1, at=l, labels=l, tick=FALSE, padj=1.1)
      }
     }
    }
        
    if(length(items)==1){
     if(add==TRUE){
      for(j in items){
       devAskNewPage(ask = TRUE)
       lines(x$monotonicity[[step]]$results[[j]][[2]][,10], type="b", ylim=ylim, xlim=c(1,dim(x$monotonicity[[N]]$results[[j]][[2]])[1]), pch=1, cex=1.5, lty=lty, lwd=lwd, col=col)
      }
     }
    }
        
    devAskNewPage(ask = FALSE)
    }

# PLOT "IRF" 
    if (type=="IRF"){
    
#
#MONOTONICITY
#
#forward plot of the item response function of the binning groups of item j -> violation of of monotonicity when IRF is decreasing
#k=4  #rest score group size
#k=10 #mean
#k=11 #P(X >=1)
#k=12 #P(X >=2)
#k=13 #P(X >=3)
#k=14 #P(X >=4)

    k <- 10
    default.ylim <- c(0,m)
    default.xlim <- c(0,N+N/40)
    #default.col <- 1
    default.lwd <- 2
    
    for(j in items){
     irf.mono <- matrix(, x$member.Rscore$n.mono[1,j], N)
     for(i in x$initial.subsample.size:N){irf.mono[,i] <- x$monotonicity[[i]]$results[[j]][[2]][,k]}
     devAskNewPage(ask = TRUE)
     plot(NA, ylim=ylim, xlim=xlim, xlab="Subsample size n", ylab=expression(italic(paste(bar(X)[j]," | ",R[(j)],sep=" "))), main=x$monotonicity[[N]]$I.labels[j], las=1)  # NEW: R[(j)] ipv R[(-j)] #
     if(n0==TRUE){axis(1, at=x$initial.subsample.size, labels=expression(italic(n[0])), tick=TRUE, padj=-1.3, cex.axis=1)}
     if(is.numeric(n1)==TRUE){axis(1, at=n1, labels=expression(italic(n[1])), tick=TRUE, padj=-1.3, cex.axis=1)}
     if(n2==TRUE){axis(1, at=x$n2, labels=expression(italic(n[2])), tick=TRUE, padj=-1.3, cex.axis=1)}
     for(bin in 1:dim(x$monotonicity[[N]]$results[[j]][[2]])[1]){
      lines(irf.mono[bin,], col=bin, lty=1, lwd=lwd)
      text(N,irf.mono[bin,N], bin, cex=1.3, pos=4)}
     }
     devAskNewPage(ask = FALSE)
    }

# PLOT "followup" 
    if (type=="followup"){
    obj <- x$objective$observation.obj
    N <- x$data$N
    default.step <- N
    default.reference.step <- min(step)-1

    quantiles <- c(.025,.25,.50,.75,.975)
    subsample.step <- x$subsample[which(x$subsample[,reference.step]>0),reference.step]
    obj.quantile <- apply(obj[subsample.step,],2,function(qu){quantile(qu,quantiles,na.rm=TRUE)})

    default.ylim <- c(0, max(obj[subsample.step,],na.rm=TRUE)+5)
    default.xlim <- c(0, N)
    
    plot(NA, xlim=xlim, ylim=ylim, ylab="Objective function", xlab="Subsample size n",las=1)
    for(i in 1:length(subsample.step)){lines(obj[subsample.step[i],], col=1, lty=i)}
    if(n0==TRUE){axis(1, at=x$initial.subsample.size, labels=expression(italic(n[0])), tick=TRUE, padj=-1.3, cex.axis=1)}
    if(is.numeric(n1)==TRUE){axis(1, at=n1, labels=expression(italic(n[1])), tick=TRUE, padj=-1.3, cex.axis=1)}
    if(n2==TRUE){axis(1, at=x$n2, labels=expression(italic(n[2])), tick=TRUE, padj=-1.3, cex.axis=1)}
    text(min(xlim),max(ylim),paste("n=",reference.step,sep=""), pos=4)
    segments(reference.step,0,reference.step,max(ylim), col=2, lwd=2)
    
    for(it in step){
     devAskNewPage(ask = TRUE)
     plot(NA,xlim=xlim,ylim=ylim, ylab="Objective function", xlab="Subsample size n", las=1)
     if(n0==TRUE){axis(1, at=x$initial.subsample.size, labels=expression(italic(n[0])), tick=TRUE, padj=-1.3, cex.axis=1)}
     if(is.numeric(n1)==TRUE){axis(1, at=n1, labels=expression(italic(n[1])), tick=TRUE, padj=-1.3, cex.axis=1)}
     if(n2==TRUE){axis(1, at=x$n2, labels=expression(italic(n[2])), tick=TRUE, padj=-1.3, cex.axis=1)}
     for(q in 1:length(quantiles)){lines(obj.quantile[q,], col=1, lty=q)}

     incl <- which(abs(sign(match(x$subsample[which(x$subsample[,it]>0),it],x$subsample[which(x$subsample[,it-1]>0),it-1],nomatch=0))-1)==1)
     for(w in 1:length(incl)){lines(obj[x$subsample[incl[w],it],], col=1, lwd=3)}
     text(min(xlim),max(ylim),paste("n = ",it,sep=""), pos=4)
     text(min(xlim),max(ylim)-(max(ylim)/10), paste("Enter:", paste(sort(x$subsample[incl,it]), sep=" ",collapse=","),sep=" ",collapse=" "), pos=4)
     segments(it, 0, it, max(ylim))
    
     if(length(incl)>1){
      excl <- sort(which(abs(sign(match((1:N)[-x$subsample[which(x$subsample[,it]>0),it]],(1:N)[-x$subsample[which(x$subsample[,it-1]>0),it-1]],nomatch=0))-1)==1))
      for(w in 1:length(excl)){lines(obj[ (1:N)[-x$subsample[,it]][excl][w] ,] ,col=1, lwd=3, lty="22")}
      text(min(xlim),max(ylim)-(2*max(ylim)/10), paste("Leave:", paste(sort((1:N)[-x$subsample[,it]][excl]), sep=" ",collapse=","),sep=" ",collapse=" "), pos=4)
     }
    }
    devAskNewPage(ask = FALSE)
    }  

# PLOT "scale" 
    if (type=="scale"){
   
    default.xlim <- c(1,N)
    default.ylim <- c(1,J)
    default.lwd <- 4
    default.lty <- 1
    default.col <- 1
       
    max.scale <- max(x$scale,na.rm=TRUE)
    scaling <- matrix(,max.scale,N)
    for(i in 1:max.scale){
     temp <- x$scale
     temp[which(temp!=i)] <- 0
     scaling[i,] <- apply(temp,2,sum)/i
    }
    
    default.scale <- 0:max.scale
    for(ids in id.scale){
    devAskNewPage(ask = TRUE)
    if(ids!=0){
     if((max.scale-ids)<0)stop(paste("Maximum number of scale is", max.scale))
     temp <- max.scale-ids+1
     interest.scale <- apply(scaling,2,order)[temp,]
     is.na(interest.scale[1:(x$initial.subsample.size-1)]) <- TRUE
    } else interest.scale <- 0
   
   plot(NA, xlim=xlim, ylim=ylim, xlab="Subsample size n",yaxt="n", ylab=paste("Items belonging to scale", ids))
   axis(2, at=items, labels=colnames(x$data$X)[items], las=1)
   if(n0==TRUE){axis(1, at=x$initial.subsample.size, labels=expression(italic(n[0])), tick=TRUE, padj=-1.3, cex.axis=1)}
   if(is.numeric(n1)==TRUE){axis(1, at=n1, labels=expression(italic(n[1])), tick=TRUE, padj=-1.3, cex.axis=1)}
   if(n2==TRUE){axis(1, at=x$n2, labels=expression(italic(n[2])), tick=TRUE, padj=-1.3, cex.axis=1)}
   for(j in items){
    temp1 <- rep(NA,N)
    temp1[which(x$scale[j,]==interest.scale)] <- j
    lines(temp1, col=col, lwd=lwd, lty=lty)
   }
   }
  devAskNewPage(ask = FALSE)
  }
 

# PLOT "num.scale" 
    if (type=="num.scale"){
      
    max.scale <- max(x$scale,na.rm=TRUE)
    scaling <- matrix(,max.scale,N)
    for(i in 1:max.scale){
     temp <- x$scale
     temp[which(temp!=i)] <- 0
     scaling[i,] <- apply(temp,2,sum)/i
    }
    tmp1 <- scaling
    is.na(tmp1[(which(tmp1==0))]) <- TRUE
    tmp1[which(is.na(tmp1)==FALSE)] <- 1
    tmp2 <- apply(tmp1,2,function(x)sum(x,na.rm=TRUE))
    is.na(tmp2[1:x$initial.subsample.size]) <- TRUE
    
    default.xlim <- c(1,N)
    default.ylim <- c(1,max.scale)
    default.lwd <- 2
    default.lty <- 1
    default.col <- 1
     
    plot(tmp2, type="l", xlim=xlim, ylim=ylim, xlab="Subsample size n", ylab="Number of scales", lab=c(5,max.scale,7), las=1, lwd=lwd, col=col, lty=lty)
    if(n0==TRUE){axis(1, at=x$initial.subsample.size, labels=expression(italic(n[0])), tick=TRUE, padj=-1.3, cex.axis=1)}
    if(is.numeric(n1)==TRUE){axis(1, at=n1, labels=expression(italic(n[1])), tick=TRUE, padj=-1.3, cex.axis=1)}
    if(n2==TRUE){axis(1, at=x$n2, labels=expression(italic(n[2])), tick=TRUE, padj=-1.3, cex.axis=1)}
   }
}
