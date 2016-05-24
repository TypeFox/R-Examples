summary.DocWordTable<-function (object, nword=50, ordFreq=TRUE, ...) 
{
  res <- object
     if (!inherits(res, "DocWordTable")) 
        stop("non convenient object") 
      cat("\nDocWordTable summary\n")
      cat("\nNumber of documents\n")
	cat(" ",res$Ndoc,"\n")
	cat("\nCorpus size\n")
	cat("",res$Nlength,"\n")
	cat("\nVocabulary size\n")
	cat("",res$Nword,"\n")
 if(ordFreq){
      nword<-min(nword,nrow(res$Tfreq))
      cat("\nGlossary of the ",nword," most frequent  words\n")
      print(res$Tfreq[c(1:nword),])
   }else{
     Tabfq <- cbind(as.data.frame(res$Nfreqword), as.data.frame(res$Ndocword))
     colnames(Tabfq) <- c("Frequency", "N.Documents")
     nword<-min(nword,nrow(Tabfq))  
     cat("\nGlossary of the ",nword," most frequent  words\n")
     print(Tabfq[c(1:nword),])   
   }
    
}
