identRow<-function(x) {
  # x: table of hypothesis (a matrix)
  # example
  #   x
  #   Black Grey White        
  #   0    1     0 
  #   0    0     1 
  #   1    0     0 
  #   1    1     1 
  
  # n: nb of rows of table x
  n<-ncol(x)
  # m : nb of lines of input matrix of hypothesis
  m<-nrow(x)
  # find non zero positions (identification of hypothesis)
  pos<-x*col(x) 
  hyp<-colnames(x)
  #  zrow<-numeric()
  if (nrow(pos) == 0) {
    return(NULL)
  } else {
    ident<-list(paste(hyp[pos[1,]],collapse=" "))
  }
  if (nrow(pos) == 1) {
    ident1<-unlist(ident)
    return(ident1)
  } else {
    for(i in 2:m)
    ident<-c(ident,paste(hyp[pos[i,]],collapse=" "))
    ident1<-unlist(ident)
    return(ident1)
  }
}
