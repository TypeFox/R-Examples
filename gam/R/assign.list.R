assign.list<-function(assignx,term.labels){
  ass<-as.list(seq(term.labels))
  names(ass)<-term.labels
  indexset<-seq(along=assignx)
  lapply(ass,function(i,indexset,assignx)indexset[assignx==i],indexset,assignx)
}
