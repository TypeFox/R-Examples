stepwise.seqs<-function(file, remGaps=TRUE) {
  info<-scan(file, n=2, quiet=TRUE)
  numseqs<-info[1]; seqlength<-info[2]
  chStrings<-scan(file,what=character(),sep="\n", quiet=TRUE) 
  #read each row of file into a character string
  chStrings<-chStrings[-1] #don't need first row
  seqNames<-rep(NA,numseqs)
  dat<-matrix(NA,nrow=numseqs,ncol=seqlength)
  for(i in 1:numseqs) {
    seqNames[i]<-substring(chStrings[i],first=1,last=10)
    dat[i,]<-substring(chStrings[i],first=(11:(10+seqlength)),
                       last=(11:(10+seqlength)))
  }
  colnames<-paste("s",(1:seqlength),sep="")
  dimnames(dat)<-list(seqNames,colnames)
  if(remGaps) {
    indmat<-dat=="-"
    chuckcol<-apply(indmat,2,any)
    dat<-dat[,!chuckcol]
}
return(dat)	
}

stepwise.seqs<-stepwise.seqs("simulinfile.txt")
save(stepwise.seqs,file="stepwise.seqs.RData")

