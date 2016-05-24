ReadSequence <-function(Seq,m,background,convertDNA){
  
  matriuSeq<-matrix(,(length(Seq)-m),m)
  for(i in 1:(length(Seq)-m)){
    matriuSeq[i,]<-Seq[i:i+m]
  }
  matriuSeq<-apply(matriuSeq,1,convertDNA=convertDNA, background=background)
  matriuSeq<-t(matriuSeq)

  matriuSeq
}

