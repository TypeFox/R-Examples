summary.SegmentsRep <-function (object, nseg=20, ordFreq=FALSE, ...) 
{
    res <- object
     if (!inherits(res, "SegmentsRep")) 
        stop("non convenient object")
	cat("\nNumber of segments \n")
	cat("",res$nbseg,"\n")
   if(ordFreq){
     cat("\nGlossary of the ",nseg," most frequent repeated segments\n")
	 nseg<-min(nseg,nrow(res$Glossary.segments$segOrderlist))
      repeated.segments<-res$Glossary.segments$segOrderFreq[1:nseg,]
   }else{
     cat("\nGlossary of the ",nseg," first repeated segments in lexicometric order\n")
	 nseg<-min(nseg,nrow(res$Glossary.segments$segOrderlist))
       repeated.segments<-res$Glossary.segments$segOrderlist[1:nseg,]
   }
       
    print(repeated.segments)
 }
