hamming.distance <-
  function(xsource, targets){
    valid.indexes<-which(nchar(targets)==nchar(xsource))
    distances<-rep(NA, length(targets))
    distances[valid.indexes]<-stringdist(xsource, targets[valid.indexes], method='hamming')
    names(distances)<-targets
    return(distances)
}
