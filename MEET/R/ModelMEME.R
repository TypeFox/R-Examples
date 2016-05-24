ModelMEME<-function(iicc){

     num_motif<-iicc$nummotif
    len_motif<-iicc$lenmotif
    factor<-as.matrix(iicc$Transcriptionfactor)
    listfactor<-lapply(c(1:nrow(factor)),function(x){factor[x,]})
    write.fasta(listfactor, names=c(1:(nrow(factor))), nbchar = ncol(factor), file.out="background.fa",open="w")
    system(paste(paste(paste(paste(iicc$meme, " background.fa -oc memeout -dna -mod anr -w", sep =" "), len_motif,sep=" "), "-nmotifs", sep=" "),num_motif,sep=" "))

}
