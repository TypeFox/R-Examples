ModeltransMEME<-function(iicc){

  call.transfac2meme<-iicc$transfac2meme
  write.fasta <- get("write.fasta",pos="package:seqinr")
  factor<-as.matrix(iicc$Transcriptionfactor)
  listfactor<-lapply(c(1:nrow(factor)),function(x){factor[x,]})
  writeMEME(factor,0)


}
