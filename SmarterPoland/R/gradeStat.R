# support function
# plots two dimensional grade analysis on one panel
plotGradeStat2D <- function(variabl1, variabl2, Xaxis = "", Yaxis = "", cex.text=0.8, addLabels=TRUE) {
  tab    <- table(factor(variabl1),factor(variabl2))
  tabSum <- addmargins(tab, 2)
  tabProp<- prop.table(tabSum, 2)
  tabCS  <- apply(tabProp, 2, cumsum)
  
  kolor  <- c("#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", 
              "#B3DE69", "#FCCDE5", "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F"
                )[1:ncol(tab)]
  plot(c(0,1),c(0,1),type="n",pch=19,xlab=Xaxis,ylab=Yaxis)
  abline(0,1,col="grey")
  abline(h=seq(0,1,0.2),col="grey95",lty=3)
  abline(v=seq(0,1,0.2),col="grey95",lty=3)
  for (i in 1:ncol(tab)) {
    points(c(0,tabCS[,"Sum"]), c(0,tabCS[,i]), type="b", pch=19, col=kolor[i])
  }
  legend("topleft", colnames(tab), col=kolor, pch=10, lwd=3,bty="n")
  
  par(xpd=NA)
  if (addLabels) 
    text(tabCS[,"Sum"], apply(tabCS,1,min),rownames(tabCS), srt=-45, adj=c(0,0),cex=cex.text, col="black")
  par(xpd=F)
}
# this function plots two panels with flipped axes
plotGradeStat2D2  <- function(variabl1, variabl2, Xaxis="", Xaxis1=Xaxis, Xaxis2=Xaxis, Yaxis="", Yaxis1=Yaxis, Yaxis2=Yaxis, ...) {
  par(mfrow=c(1,2))
  par(xpd=F)
  plotGradeStat2D(variabl1, variabl2, Xaxis=Xaxis1, Yaxis=Yaxis1, ...)
  plotGradeStat2D(variabl2, variabl1, Xaxis=Xaxis2, Yaxis=Yaxis2, ...)
}

plotGradeStat <- function(variabl1, variabl2, decreasing = TRUE, Xaxis = "", Yaxis = "", skala=c(0.005,0.5), cex.text=0.8, cutoff = 0.01) {
  zm1r   <- variabl1/sum(variabl1)
  zm2r   <- variabl2/sum(variabl2)
  iloraz <- zm1r/zm2r
  if (decreasing) {
    zm1r   <- zm1r[order(iloraz, decreasing=FALSE), 1, drop=FALSE]
    zm2r   <- zm2r[order(iloraz, decreasing=FALSE), 1, drop=FALSE]
    iloraz <- zm1r/zm2r
  }
  par(mfrow=c(1,2))
  par(xpd=F)
  # first plot is the grade analysis
  # for one dimensional data
  plot(c(0,cumsum(zm1r[,1])),c(0,cumsum(zm2r[,1])),type="b",pch=19,xlab=Xaxis,ylab=Yaxis)
  abline(0,1,col="grey")
  par(xpd=NA)
  # move labels apart
  odleglosci <- sqrt(diff(c(0,cumsum(zm1r[,1])))^2+diff(c(0,cumsum(zm2r[,1])))^2)
  korekta    <- numeric(length(odleglosci))
  for (i in seq_along(korekta)) {
    if (odleglosci[i] < cutoff) 
      korekta[i] <- cutoff + korekta[i-1]
  }
  text(cumsum(zm1r[,1])+korekta+2*cutoff,cumsum(zm2r[,1])+korekta-2*cutoff,rownames(zm1r), srt=-45, adj=c(0,0),cex=cex.text)
  text(cumsum(zm1r[,1])+korekta-2*cutoff,cumsum(zm2r[,1])+korekta+2*cutoff,paste(round((1/iloraz[,1]-1)*1000)/10," %",sep=""), srt=-45, adj=c(1,1),cex=cex.text)
  par(xpd=F)

  plot(1,type="n",log="xy",xlim=skala,ylim=skala, las=1, cex.axis=0.8, xlab=Xaxis, ylab=Yaxis)
  abline(0,1,col="grey")
  abline(h=c(0.0001*c(1,2,5),0.001*c(1,2,5),0.01*c(1,2,5),0.1*c(1,2,5)),col="grey95")
  abline(v=c(0.0001*c(1,2,5),0.001*c(1,2,5),0.01*c(1,2,5),0.1*c(1,2,5)),col="grey95")
  points(zm1r[,1],zm2r[,1],pch=19)
  par(xpd=NA)
  text(zm1r[,1],zm2r[,1],rownames(zm1r), srt=-45, adj=c(-0.1,-0.1),cex=cex.text)
  par(xpd=F)
}

