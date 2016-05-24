getConstraint_rec2 <-
function(seqsf1,frame,code,n,cpt,pepList){

  ## In overlapping frame (f2)
  
  # If anti-sense overlap, generate seqsf2 the reverse complement of dna sequences seqsf1.
  if(frame <1){
    seqsf2=compRev(as.matrix(seqsf1))
  }else{
    seqsf2=seqsf1
  }
  
  # Remove peptides with stops (impossible in double coding sequences)
  pep2=deleteStops(unique(toPep(as.array(seqsf2),posLecturef2(frame))))
  
  # Get the number assoiated with the overlapping peptides (position of the peptides in pepList[[n]] + 20^n )
  l2=which(pepList[[n]] %in% pep2)+20^n

  # If l2 is not empty, generate allSeqsf2, the set of the dna sequences encoding for the overapping protein in the overlapping frame
  if(length(l2)>0){
    allSeqsf2=unlist(apply(as.matrix(pep2),1,getAllSeqs,posLecturef2(frame),code))

    # Remove dna sequences already considered (directly from seqsf1)
    allSeqsf2=allSeqsf2[which(!(allSeqsf2 %in% seqsf2))] 

    if(length(allSeqsf2)>0){
      cpt=cpt+1
      # Recursive call to getConstraint_rec2 for the graph traversal
      l2=c(l2,getConstraint_rec1(allSeqsf2,frame,code,n,cpt,pepList))
    }
  }
 
  return(l2)
}
