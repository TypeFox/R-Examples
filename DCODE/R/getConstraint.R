getConstraint <-
function(pep,frame,code,pepList){
  
  n=nchar(pep)
  
  if(n>5){
    print("WARNING : Please enter a peptide with length less than or equal to 5 amino acids.")
  }else{
    
  cpt=0
  l=which(pepList[[n]] == pep)
  seqsf1=getAllSeqs(pep,posLecturef1(frame),code)

  # Recursive call to getConstraint_rec2 for the graph traversal
  res = sort(c(l,getConstraint_rec2(seqsf1,frame,code,n,cpt,pepList)))

  print("CAPITAL LETTERS and small letters denote respectively peptides read in the REFERENCE frame and the overlapping frame.")
  
  # Print the associated constraint
  printConstraint(res,n,pepList)
  }
}
