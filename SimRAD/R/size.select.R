size.select <-
function(sequences, min.size, max.size, graph=TRUE, verbose=TRUE){
  ssel <- sequences[width(sequences) < max.size & width(sequences) > min.size]
  if(verbose==TRUE){
    cat(length(ssel), " fragments between ", min.size, " and ", max.size, " bp \n", sep="")
  }
  if(graph==TRUE){
    bk<-hist(width(sequences), breaks=length(sequences)/20, plot=FALSE)$breaks
    hist(width(sequences),border="grey75", col="grey75", breaks=bk, main="", xlab="Locus size (bp)", ylab="Number of loci")
    hist(width(ssel),border="red", col="red", add=TRUE, breaks=bk)
    text(mean(c(min.size, max.size)), max(hist(width(ssel), breaks=bk, plot=FALSE)$counts), pos=4, labels=paste(length(ssel), " loci between ", min.size, " and ", max.size, " bp", sep=""), col="red", cex=0.9, font=2)
  }
  return(ssel)
}
