# -------------------------------------------------------------------------------------
# this function performs peptide to protein roll-up
# -------------------------------------------------------------------------------------

# arguments ___________________________________________________________________________
#           : pep.Expr.Data   - matrix of peptide expression data
#           : rollup.map      - the map to peptide to protein mapping

# output ______________________________________________________________________________
#           : prot.Expr.Data  - matrix of peptide expression data 

pep2prot = function(pep.Expr.Data, rollup.map){
  
  protIDs = unique(rollup.map)
  prot.Expr.Data = matrix(,nrow = length(protIDs),ncol = dim(pep.Expr.Data)[2])
  
  for (i in 1:length(protIDs)){
    
    temp = protIDs[i]
    pepSet = pep.Expr.Data[which(rollup.map==temp),]
    
    if (length(which(rollup.map==temp))==1){
      prot.Expr.Data[i,] = pepSet
    }
    else{
      prot.Expr.Data[i,] = apply(pepSet,2,median,na.rm = T)
    }
    
  }
    
  return(prot.Expr.Data)
}