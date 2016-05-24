multtest.cor <- function(mult.var,uni.var,method="pearson",p.method="fdr",ordered=TRUE) {
  nvar <- ncol(mult.var)
  tab.res <- data.frame(Corr=integer(nvar),P.value=integer(nvar)," "=character(nvar),
    row.names=colnames(mult.var),check.names=FALSE,stringsAsFactors=FALSE)
  for (i in 1:nvar) {
    test <- suppressWarnings(cor.test(uni.var,mult.var[,i],method=method))
    tab.res$Corr[i] <- test$estimate
    tab.res$P.value[i] <- test$p.value
  }
  tab.res$P.value <- p.adjust(tab.res$P.value,method=p.method)
  tab.res$Corr <- signif(tab.res$Corr,4)
  tab.res$P.value <- signif(tab.res$P.value,5)
  tab.res[,3] <- .psignif(tab.res$P.value)
  if (ordered) {tab.res <- tab.res[order(tab.res$Corr,decreasing=TRUE),]}
  res <- list(tab=tab.res,p.method=p.method)
  class(res) <- c("multtest","multtest.cor","list")
  return(res)
}

plot.multtest.cor <- function(x,arrows=TRUE,main=NULL,pch=16,cex=1,
  col=c("red","orange","black"),labels=NULL,...) {
  to.print <- as.data.frame(x$tab$Corr)
  rownames(to.print) <- rownames(x$tab)
  if (length(col)<3) {col <- rep(col[1],3)}
  cols <- rep(col[3],nrow(x$tab))
  cols[which(x$tab$P.value<0.05)] <- col[1]
  cols[which(x$tab$P.value>0.05 & x$tab$P.value<0.1)] <- col[2]
  MVA.corplot.1comp(to.print,fac=NULL,"Correlation",
    arrows=arrows,main=main,pch=pch,cex=cex,col=cols,lwd=1,labels=labels,legend.col=NULL,
    drawintaxes=TRUE)
}



