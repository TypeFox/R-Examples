split.direct.sum<-function(x){  
  
  if(is.na(x)){
   return(NULL)
  }else{
    x <- gsub("\n", "", x, fixed=TRUE)  # required if the random formula is very long and broken over lines
    openB<-gregexpr("\\(", x)[[1]]
    closeB<-gregexpr("\\)", x)[[1]]
    true_openB<-openB
    true_closeB<-closeB

    for(i in 1:length(openB)){
      dist<-outer(openB, closeB, function(x,y){y-x})
      dist[which(dist<0)]<-Inf
      dist<-which(dist==min(dist), arr.ind=T)[1,] 
      true_openB[i]<-openB[dist[1]]
      true_closeB[i]<-closeB[dist[2]]
      openB<-openB[-dist[1]]
      closeB<-closeB[-dist[2]]
    }

    plus<-gregexpr("\\+", x)[[1]]
    internals<-matrix(mapply(function(x,y){(plus>x & plus<y)}, x=true_openB, y=true_closeB), length(plus), length(true_openB))
    rterms<-strsplit(x, "")[[1]]
    rterms[plus[which(rowSums(internals)!=0)]]<-"leaveMCMCleave"
    rterms<-paste(rterms, collapse="")
    rterms<-strsplit(rterms, " *\\+ *")[[1]]
    rterms<-gsub("leaveMCMCleave", "+", rterms)
    return(rterms)
  }
}

