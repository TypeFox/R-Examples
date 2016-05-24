 mruns.out <- function(prefix, lab.out="", input, param, output, tXlab, pdfout=F){

 # produces output files for multiple runs and one variable

 # produce graph 
  t <- output[,1]; X <- output[,-1]
  if(pdfout==T) pdf(file=paste(prefix,"-",lab.out,"-out.pdf",sep=""))
   par(xaxs="i", yaxs="i")
   matplot(t, X, type="l", col=1, xlab=tXlab[1], ylab=tXlab[2])
   legend ("top", legend=paste(param$plab,"=",param$pval),
           lty=1:length(param$pval),col=1,cex=0.8)
  if(pdfout==T) dev.off()

 
 # writes output file as csv
  fileout.csv=paste(prefix,"-",lab.out,"-out.csv",sep="")
  write.table(input,fileout.csv, sep=",", row.names=F, quote=F)
  name.param <- t(c(param$plab,param$pval)) 
  write.table(name.param, fileout.csv, sep=",",
              col.names=F, row.names=F, append=T, quote=F)
  name.out <- t(names(output))
  write.table(name.out,fileout.csv, sep=",", col.names=F, 
              row.names=F, append=T, quote=F)
  write.table(output,fileout.csv, sep=",", col.names=F, row.names=F, 
              append=T, quote=F)

}
