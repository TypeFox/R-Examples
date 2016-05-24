PlotFDRs <- function(lpcfdr.out, frac=.25){
  # frac is the fraction of genes plotted.
  CheckPlotFDRsFormat(lpcfdr.out,frac)
  lpcfdr <- lpcfdr.out$fdrlpc
  tfdr <- lpcfdr.out$fdrt
  tfdrs <- tfdr[tfdr<quantile(tfdr,frac)]
  lpcfdrs <- lpcfdr[lpcfdr<quantile(lpcfdr,frac)]
  plot(tfdrs[order(tfdrs, decreasing=FALSE)], type="l", ylim=range(c(tfdrs,lpcfdrs)), main="Estimated FDRs for T and LPC", xlab="Num Genes Called Significant", ylab="Estimated FDR")
  points(lpcfdrs[order(lpcfdrs,decreasing=FALSE)], type="l", col="red")
  legend("topleft", pch=15, col=c("black", "red"), c("T", "LPC"))
}
  
