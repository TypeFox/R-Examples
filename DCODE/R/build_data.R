build_data <-
function(){
  print("Initialization...")
  
  # List of the amino acids
  aaList<-sort(unique(SEQINR.UTIL$CODON.AA[,3]))[-1]
  
  # Standard genetic code
  tmp<-list()
  for(i in aaList){
    tmp[[i]]<-as.matrix(SEQINR.UTIL$CODON.AA)[which(SEQINR.UTIL$CODON.AA[,3]==i),1]
  }
  geneticCode<-tmp
 

  # Build peptideList
  diPep<-words(length=2, alphabet=aaList)
  triPep<-words(length=3, alphabet=aaList)
  quadriPep<-words(length=4, alphabet=aaList)
  pentaPep<-words(length=5, alphabet=aaList)
  peptideList<-list(aa=aaList,diPep=diPep,triPep=triPep,quadriPep=quadriPep,pentaPep=pentaPep)
  
  return(list(aaList=aaList, geneticCode=geneticCode, peptideList=peptideList))
}
