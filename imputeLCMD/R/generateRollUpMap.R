# -------------------------------------------------------------------------------------
# this function generates a map for peptide to protein roll-ul
# -------------------------------------------------------------------------------------

# arguments ___________________________________________________________________________
#           : nProt           - number of proteins to map to the peptide expression data
#           : pep.Expr.Data   - matrix of peptide expression data

# output ______________________________________________________________________________
#           : pep.prot.Map    - the peptide to protein map (for each row in 
#                               pep.prot.Map the corresponding value corresponds 
#                               to the index of the protein that peptide is mapped to)

  generate.RollUpMap = function(nProt, pep.Expr.Data){
    
    n = dim(pep.Expr.Data)[1]
    temp = 1:nProt
    pep.prot.Map = rep(0,n)
    pep.prot.Map[sample(temp)] = sample(temp)
    pep.prot.Map[which(pep.prot.Map==0)] = sample.int(nProt, 
                                                      size = (n - nProt), 
                                                      replace = T)
    
    return(pep.prot.Map)
  }