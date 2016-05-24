triangle.plot <-
function(hi.index=NULL, int.het=NULL, pdf=TRUE, out.file="tri_plot.pdf"){
  if (is.null(hi.index)==TRUE | is.null(int.het)==TRUE)
    stop("error, input data are not provided")
  if (is.data.frame(hi.index)==TRUE)
    hi.index<-hi.index[,2]
  if (pdf==TRUE)
    pdf(file=paste(out.file))
  plot(hi.index,int.het,xlab="Hybrid index",ylab="Interspecific heterozygosity",xlim=c(0,1),ylim=c(0,1))
  lines(c(0,0.5),c(0,1))
  lines(c(0.5,1),c(1,0))
  if (pdf==TRUE)
    dev.off()
}

