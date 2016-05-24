sample.ds<-function(dat,dur){
  nrow <- nrow(dat)
  l <- length(table(dat$dur))
  
  sample<-sample(1:nrow,nrow,replace=TRUE)
  
  l.sample<-length(table(dat[sample,]$dur))
  while((l.sample!=l)){
    
    sample<-sample(1:nrow,nrow,replace=TRUE)
    l.sample<-length(table(dat[sample,]$dur))

    }

  return(dat[sample,])
  }