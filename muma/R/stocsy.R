stocsy <-
function(scaling, threshold=FALSE, pos.threshold, neg.threshold) {
  pwd.n = paste(getwd(), "/Preprocessing_Data_", scaling, "/ProcessedTable.csv", sep="")
  x <- read.csv(pwd.n, sep=",", header=TRUE)
  x.x <- x[,2:ncol(x)]
  rownames(x.x) <- x[,1]
  x.t <- x.x
  mycor = cor(x.t, method=c("pearson"))
  library(gplots)
  col = colorpanel(50, "blue", "white", "red")
  image(mycor, axes=FALSE, col=col, main="STOCSY")
  axis(side=1, labels=colnames(mycor), at=seq(0,1,length=length(colnames(mycor))), las=2, cex.axis=0.7)
  axis(side=2, labels=colnames(mycor), at=seq(0,1,length=length(colnames(mycor))), las=2, cex.axis=0.7)
  dirout = paste(getwd(), "/STOCSY/", sep="")
  dir.create(dirout)
  o = paste(dirout, "STOCSY.pdf", sep="")
  dev.copy2pdf(file=o)
  o.cor = paste(dirout, "CorrelationMatrix.csv", sep="")
  write.csv(mycor, file=o.cor)
  if (threshold) {
    dev.new()
    image(mycor, axes=FALSE, zlim=c(pos.threshold,1), col="red", main = paste("STOCSY <", neg.threshold, " & >", pos.threshold, sep=""))
    image(mycor, axes=FALSE, zlim=c(-1,neg.threshold), col="navy", add=TRUE)
    axis(side=1, labels=colnames(mycor), at=seq(0,1,length=length(colnames(mycor))), las=2, cex.axis=0.7)
    axis(side=2, labels=colnames(mycor), at=seq(0,1,length=length(colnames(mycor))), las=2, cex.axis=0.7)
    out = paste(dirout, "STOCSY_", pos.threshold, "_", neg.threshold, ".pdf", sep="")
    dev.copy2pdf(file=out)
  }
}
