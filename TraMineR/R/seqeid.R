## ========================================
## Return a id of given subsequences
## ========================================

seqeid<-function(s){
  tmrsequenceid.internal<-function(s){
    if(is.seqe(s)){
      return(.Call(TMR_tmrsequencegetid, s))
    }
    return(NA)
  }
  #message("Event sequence analysis module is still experimental")
  if(is.seqelist(s)){
    return(as.integer(sapply(unlist(s),tmrsequenceid.internal)))
  }else if(is.seqe(s)){
    return(tmrsequenceid.internal(s))
  }else{
    stop("seq should be a seqelist. See help on seqecreate.")
  }
  return(NA)
}