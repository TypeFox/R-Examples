seequestions<-function(QB)
  {

    uq = unlist(QB)
    w = which(names(uq)=="Q")
    ques = as.vector(uq[w]) 
   # cat(ques, sep="\n")
    print(ques)
  }
