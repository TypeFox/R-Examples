hasSeg <-
function(CPvector, start, end){
  if(!(start %in% CPvector & end %in% CPvector)) return(FALSE)
  else {
    if((end != start+1) & (sum((start+1):(end-1) %in% CPvector)!=0)) return(FALSE)
    else return(TRUE)
  }
}
