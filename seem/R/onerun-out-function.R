 onerun.out <- function(prefix, lab.out="", input, output, pdfout=F){

 # produces output files for one run and one variable

 # plot the graph; in pdf output file if desired
  t <- output[,1]; X <- output[,-1]
  tlab <- names(output)[1]; Xlab <- names(output)[2]
  if(pdfout==T) pdf(file=paste(prefix,"-",lab.out,"-out.pdf",sep=""))
   par(xaxs="i",yaxs="i")
   plot(t, X, type="l", ylim=c(0,max(X)),xlab=tlab,ylab=Xlab)
  if(pdfout==T) dev.off()

 # writes output file as csv
  fileout.csv <- paste(prefix,"-",lab.out,"-out.csv",sep="")
  write.table(input,fileout.csv, sep=",", row.names=F, quote=F)
  name.out <- t(names(output))
  write.table(name.out,fileout.csv, sep=",", col.names=F, 
              row.names=F, append=T, quote=F)
  write.table(output,fileout.csv, sep=",", row.names=F,
              col.names=F, append=T,quote=F)
} 
