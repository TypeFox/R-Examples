add.value.labels<-function(x,value.labels) {
 attr(x,"value.labels")<-sort(unique(x))
 lenvallab<-length(attr(x,"value.labels"))
 if(length(value.labels) > lenvallab) {
  cat("More value labels than values, only the first",lenvallab,"will be used\n")
  value.labels<-value.labels[1:lenvallab]
 }
 names(attr(x,"value.labels"))<-value.labels
 return(x)
}
