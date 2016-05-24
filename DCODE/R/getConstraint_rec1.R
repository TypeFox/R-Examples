getConstraint_rec1 <-
function(seqsf2,frame,code,n,cpt,pepList){
  ## In reference frame (f1)
  
  # If anti-sense overlap, generate seqsf1 the reverse complement of dna sequences seqsf2.
  if(frame <1){
    seqsf1=compRev(as.matrix(seqsf2))
  }else{
    seqsf1=seqsf2
  }
  
  # Remove peptides with stops (impossible in double coding sequences)
  pep1=deleteStops(unique(toPep(as.array(seqsf1),posLecturef1(frame))))
  
  # Get the position of the peptides in pepList[[n]]
  l1=which(pepList[[n]]%in%pep1)
  
    # If l1 is not empty, generate allSeqsf1, the set of the dna sequences encoding for the overapping protein in the reference frame
    if(length(l1)>0){
    allSeqsf1=unlist(apply(as.matrix(pep1),1,getAllSeqs,posLecturef1(frame),code))

    # Remove dna sequences already considered (directly from seqsf2)
    allSeqsf1=allSeqsf1[which(!(allSeqsf1 %in% seqsf1))] ## delete sequences already seen in seqsf2

    if(length(allSeqsf1)>0){
      cpt=cpt+1
      # Recursive call to getConstraint_rec2 for the graph traversal
      l1=c(l1,getConstraint_rec2(allSeqsf1,frame,code,n,cpt,pepList))
    }
  }

  return(l1)
}
