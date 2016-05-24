

copy.attributes <- function(from,to,delete=c('names','row.names','class','dim','dimnames'),delete2=NULL) {
  a <- attributes(from)
  a[c(delete,delete2)] <- NULL
  attributes(to) <- c(attributes(to),a)
  return(to)
}

"copy.attributes<-" <- function(to,delete=c('names','row.names','class','dim','dimnames'),delete2=NULL,value) {
  return(copy.attributes(from=value,to=to,delete=delete,delete2=delete2))
}