deleteStops <-
function(vec){
  ## if peptide length >1
  if(nchar(vec)[1]>1){
    if(length(vec)==1){
      if(sum(s2c(vec)=="*")>0){vec=vec[-1]}
      return(vec)
    }
    return(vec[which(apply(is.na(apply(apply(as.array(vec),1,s2c),2, match,"*")),2,sum)==nchar(vec[1]))])    
  }else{
    ## if peptide length = 1
    return(vec[vec!="*"])
  }
}
