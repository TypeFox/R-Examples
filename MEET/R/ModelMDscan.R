ModelMDscan<-function(iicc){
  write.fasta <- get("write.fasta",pos="package:seqinr")
  read.fasta <- get("read.fasta",pos="package:seqinr")
  call.MDscan<-iicc$MDscan
  len_motif<-iicc$lenmotif
  num_motif<-iicc$nummotif
  factor<-as.matrix(iicc$Transcriptionfactor)
  listfactor<-lapply(c(1:nrow(factor)),function(x){factor[x,]})
  write.fasta(input=listfactor,names="sequenciaEstudi", nbchar = ncol(factor), file.out="sequenciaEstudi.fa",open="w")
  mdscan<-paste(call.MDscan," -i sequenciaEstudi.fa -w")
  system(paste(paste(paste(paste(paste(paste(paste(mdscan, len_motif,sep=" "),"-t", sep=" "), num_motif,sep=" "),"-b background.fa", sep=" "),"-r",sep=" "), num_motif,sep=" ") ,"-o res")) 
    
}
