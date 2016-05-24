#' volcano plot
#' @param foldchange - fold change values
#' @param pvals pvalues
#' @param pthresh pvalue threshold
#' @param ratiothresh threshold of foldchange
#' @param labels - optional labels
#' @param cex size of labels
#' @param xlab - x axis label
#' @param ylab - y axis label
#' @param xlim - xlim
#' @param cex.point - point size
#' @param main - main title
#' @export
#' @examples
#' foldchange <- rnorm(1000)
#' pval <-rexp(1000)
#' abline(v=0.05,col=2)
#'
#' volcanoplot(foldchange, pval,pthresh=0.1, ratiothresh=1,main='test')
#'
volcanoplot = function(foldchange,
                       pvals ,
                       pthresh = 0.05,
                       ratiothresh = 2,
                       xlab ="log2(T/N)" ,
                       ylab = "-log10(P)",
                       labels = NULL,
                       cex=0.6,
                       cex.point=1,
                       xlim=NULL,
                       main=NULL
){
  dataframe <- data.frame("ratio" = foldchange, "pvals" = pvals )
  rownames(dataframe) = labels

  bla = tryCatch( plot(dataframe$ratio,-log10(dataframe$pvals),col="#00000033",pch=19,xlab=xlab, ylab=ylab,xlim=xlim,main=main),
                  warning=function(bla){ dev.off(); return(1)}
  )
  if(!is.null(bla)){
    plot(dataframe$ratio,-log10(dataframe$pvals),col=1,pch=19,xlab=xlab, ylab=ylab ,xlim=xlim,cex=cex.point,main=main)
  }


  upsubset<-subset(dataframe,pvals < pthresh & ratio > ratiothresh)
  points(upsubset$ratio,-log10(upsubset$pvals),col=2,pch=19)
  points(upsubset$ratio,-log10(upsubset$pvals),col=1,pch=1)
  if(length(rownames(upsubset)) > 0){
    text(upsubset$ratio, -log10(upsubset$pvals),rownames(upsubset),cex=cex,pos=4)
  }

  abline(h=-log10(pthresh),col="gray")
  downsubset<-subset(dataframe,pvals<pthresh & ratio < -ratiothresh)
  points(downsubset$ratio,-log10(downsubset$pvals),col=3,pch=19)
  points(downsubset$ratio,-log10(downsubset$pvals),col=1,pch=1)
  if(length(rownames(downsubset)) > 0){
    text(downsubset$ratio, -log10(downsubset$pvals),rownames(downsubset),cex=cex,pos=2)
  }

  abline(v=c(-ratiothresh,ratiothresh),lty=2)
  abline(v =0,lty=2,lwd=1.5)
  return(list(upsubset=upsubset,downsubset=downsubset))
}
