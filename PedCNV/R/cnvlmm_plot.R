##' Makes formatted plots from the clustering result returned from \code{\link{ClusProc}}.
##'
##' @title Plots clustering result
##' @param x The clustering results obtained from \code{\link{ClusProc}}.
##' @param type Factor. For specifying the plot type. It must be one of 'histo', 'scat' and 'sil'. If it is 'histo', the histogram is obtained with the first PC score of the intensity measurement. For 'scat', the first PC score of the intensity measurement is plotted against the mean of the intensity measurement. For 'sil', the silhouette score is plotted. See details.
##' @param adjust Logicals. If TRUE (default), the silhouette-adjusted clustering result will be used. If FALSE, the initial clustering result will be used. See details in \code{\link{ClusProc}}.
##' @param ... Usual arguments passed to the qplot function.
##' @details
##' \itemize{
##' \item{type}{We provide three types of plots: 'hist', 'scat' and 'sil'. The first two plots are used to visually check the performance of clustering. Different clusters are represented by using different colors. The 'sil' plot is the the overview of the silhouette value for all the individuals, the silhouettes of the different clusters are printed below each other. The higher silhouettes value means the better performance.}
##' }
##' @author Meiling Liu
##' @method plot clust
##' @examples
##' # Fit the data under the given clustering numbers
##' clus.fit <- ClusProc(signal=signal,N=2:6,varSelection='PC.9')
##' plot(clus.fit,type='histo')
##' @export
plot.clust <- function(x,type=c('histo','scat','sil'), adjust=TRUE, ...){

    MEAN=NULL
    PC1=NULL
    type <- match.arg(type)
    if(type=='histo'){
        sX <- as.matrix(x$signal)
        pca <- princomp(sX)
        PCA1 <- sX%*%loadings(pca)[,1]
        if(adjust) {
            clusters <- matrix(x$silWidth$adjusted$silRes.adjust$clus)
            rownames(clusters) <- rownames(x$sil$adjusted$silRes.adjust)} else {
                clusters <- matrix(x$silWidth$unadjusted$silRes$clus)
                rownames(clusters) <- rownames(x$sil$unadjusted$silRes)}
        temp <- as.data.frame(merge(PCA1,clusters,by='row.names')[,-1])
        colnames(temp) <- c('PC1','clusters')
        temp[,2] <- factor(temp[,2])
        print(qplot(PC1,fill=clusters,data=temp,geom="density")+geom_histogram(aes_string(y='..count..'),binwidth=0.2))
  }
    
  if(type=='scat'){

      signal <- x$signal
      sX <- as.matrix((signal))
      pca <- princomp(sX)
      PCA1 <- sX%*%loadings(pca)[,1]
      segmean <- as.matrix(apply(sX,1,mean))

      if(adjust) {
          clusters <- matrix(x$silWidth$adjusted$silRes.adjust$clus)
          rownames(clusters) <- rownames(x$sil$adjusted$silRes.adjust)} else {
              clusters <- matrix(x$silWidth$unadjusted$silRes$clus)
              rownames(clusters) <- rownames(x$sil$unadjusted$silRes)}
      

      temp <- as.data.frame(merge(cbind(segmean,PCA1),clusters,by='row.names')[,-1])
      colnames(temp) <- c('MEAN','PC1','clusters')
      temp[,3] <- factor(temp[,3])
      print(qplot(MEAN,PC1,color=clusters,data=temp))
  }

  if(type=='sil'){

      if(adjust){
          silRes <- x$silWidth$adjusted$silRes.adjust
          silMean <- x$silWidth$adjusted$silMean.adjust
          clusNo <- x$silWidth$adjusted$clusNum.adjust
          clusAvg <- x$silWidth$adjusted$clusAvg.adjust 
          abandon_num <- length(x$silWidth$adjusted$abandon.id)
      }else{
          silRes <- x$silWidth$unadjusted$silRes
          silMean <- x$silWidth$unadjusted$silMean
          clusNo <- x$silWidth$unadjusted$clusNum
          clusAvg <- x$silWidth$unadjusted$clusAvg
          abandon_num <- 0
      }
      silRes <- silRes[with(silRes,order(silRes$clus,-silRes$sil)),]
      obsNo <- dim(silRes)[1]
      clusRes <- silRes$clus
      s <- rev(silRes[,"sil"])
      space <- c(0,rev(diff(cli <- silRes$clusRes)))
      space[space!=0] <- 5
      xlab <- expression("Silhouette width"* s[i])
      main <- paste("Silhouette plot")
      sub <- paste("Average silhouette width:",round(silMean,4))
      y <- barplot(s,width=1,space=space,xlim=c(min(0,min(s)),1),horiz=TRUE,col="grey",mgp=c(2.5,1,0),las=1,border=0,xlab=xlab)
      title(main=main,sub=sub,adj=0)
      mtext(paste("n=",obsNo,'; abandon=', abandon_num),adj=0)
      mtext(substitute(k ~ ~"clusters" ~ ~C[j], list(k=clusNo)),adj=1)
      mtext(expression(paste(j, " :  ", n[j], " | ", ave[i %in% Cj] ~ ~s[i])), adj = 1.04, line = -1.2)
      y <- rbind(rev(y),(clusRes))
      for (j in 1:clusNo) {
          yj <- mean(y[1,y[2,]==j-1])
          text(1,yj , paste(j-1, ":  ",table(clusRes)[j], " | ", format(clusAvg[j], digits = 1, nsmall = 2)), xpd = NA, adj = 0.8)
      }
  }
}
