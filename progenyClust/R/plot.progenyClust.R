plot.progenyClust <-
function(x,data=NULL,k=NULL,errorbar=FALSE,xlab='',ylab='',...){
  if(is.null(data)){
    if(x$method=='gap'){
      if(errorbar==T & dim(x$score)[1]>1){
        Hmisc::errbar(x$ncluster[-c(1,length(x$ncluster))],x$mean.gap,yplus=x$sd.gap+x$mean.gap,yminus=x$mean.gap-x$sd.gap,xlab=xlab,ylab=ylab,...)
      }else{
        plot(x$ncluster[-c(1,length(x$ncluster))],x$mean.gap,xlab=xlab,ylab=ylab,...)
      }
    }else if(x$method=='score'){
      if(errorbar==T & dim(x$score)[1]>1){
        Hmisc::errbar(x$ncluster,x$mean.score,yplus=x$sd.score+x$mean.score,yminus=x$mean.score-x$sd.score,xlab=xlab,ylab=ylab,...)
      }else{
        plot(x$ncluster,x$mean.score,xlab=xlab,ylab=ylab,...)
      }
    }else{
      par(mfrow=c(2,1))
      if(errorbar==T & dim(x$score)[1]>1){
        Hmisc::errbar(x$ncluster[-c(1,length(x$ncluster))],x$mean.gap,yplus=x$sd.gap+x$mean.gap,yminus=x$mean.gap-x$sd.gap,xlab=xlab,ylab='gap',xlim=c(min(x$ncluster),max(x$ncluster)),...)
        Hmisc::errbar(x$ncluster,x$mean.score,yplus=x$sd.score+x$mean.score,yminus=x$mean.score-x$sd.score,xlab=xlab,ylab='score',...)
      }else{
        plot(x$ncluster[-c(1,length(x$ncluster))],x$mean.gap,xlab=xlab,ylab='gap',xlim=c(min(x$ncluster),max(x$ncluster)),...)
        plot(x$ncluster,x$mean.score,xlab=xlab,ylab='score',...)
      }
    }
  }else{
    if(is.null(k)){
      if(!is.null(x$mean.score)){
        if(x$score.invert==T){
          k=which(x$mean.score==min(x$mean.score))
        }else{
          k=which(x$mean.score==max(x$mean.score))
        }
      }else{
        if(x$score.invert==T){
          k=which(x$mean.gap==min(x$mean.gap))+1
        }else{
          k=which(x$mean.gap==max(x$mean.gap))+1
        }
      }
    }else{
      if(!k%in%x$ncluster){
        stop('cluster number not valid')
      }else{
        k=which(x$ncluster==k)
      }
    }
    if(nrow(x$cluster)!=nrow(data)){
      stop("size of clustering result does not fit size of dataset")
    }
    if(ncol(data)<2){
      stop("cannot plot 1D data set")
    }
    if (ncol(data) == 2) {
      xlim <- c(min(data[, 1]), max(data[, 1]))
      ylim <- c(min(data[, 2]), max(data[, 2]))
      plot(x =NULL, y = NULL, xlim = xlim, ylim = ylim, 
           xlab = xlab, ylab = ylab, ...)
      cols <- rainbow(max(x$cluster[,k]))[x$cluster[,k]]
      points(data, col = cols, pch = 19, cex = 0.8)
    }else{
      if (ncol(data) > 20) 
        stop("cannot plot more than 20 features at once")
      cols <- rainbow(max(x$cluster[,k]))[x$cluster[,k]]
      clustPanel <- function(x, y, ...) {
        points(x, y, col = cols, pch = 19, cex = 0.8)
      }
      if(length(colnames(data))>0){
        labels=colnames(data)
      }else{
        labels=paste0(deparse(substitute(data, env = parent.frame())),"[,",1:ncol(data),"]",sep="")
      }
      pairs(data, labels, lower.panel = clustPanel, upper.panel = clustPanel, ...)
    }
  }
  
}
