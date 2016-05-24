 rnum.out <- function(prefix, lab.out="", input, output, tXlab, pdfout=F){

 # produces output files for n realization of one run and one variable


 # plot the graph; in pdf output file if desired
  t <- output[,1]; X <- output[,-1]
 # Upper and lower bounds 
  Xm <- X[,1]; Xs <- X[,2]
  
  tlab <- tXlab[1]; Xlab <- tXlab[2]

  if(pdfout==T) pdf(file=paste(prefix,"-",lab.out,"-out.pdf",sep=""))
   par(xaxs="i",yaxs="i")
   plot(t, Xm, type="l", ylim=c(0,max(X+Xs)),xlab=tlab,ylab=Xlab)
   lines(t, Xm-Xs, lty=2)
   lines(t, Xm+Xs, lty=4)
   legend("top",legend=c("Avg", "Avg-sd", "Avg+sd"),lty=c(1,2,4),col=1)
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
