


mobility_index = function(seq,alphabet,nalpha){
  cumulative=1
  
  for(j in 1:(length(seq)-1)) {
    if(seq[j]!=seq[j+1]) {
      cumulative=cumulative+1
    }
  }
  
  return(cumulative)
}


