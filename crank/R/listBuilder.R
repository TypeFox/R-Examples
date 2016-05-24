# Builds a list by descending x until it finds a non-list element,
# then applies FUN to that element if FUN is not NULL
# Return value:
# If FUN was not invoked, the list structure produced by the
# recursive calls to listBuilder.
# If FUN has been invoked, its return value, which becomes an
# element of the current list structure

listBuilder<-function(x,FUN=NULL,fargs=NULL) {
 
 if(is.list(x)) {
  lenx<-length(x)
  newlist<-vector("list",lenx)
  for(lindex in 1:lenx) {
   newlist[[lindex]]<-listBuilder(x[[lindex]],FUN,fargs)
  }
  return(newlist)
 }
 else {
  # if FUN is NULL, it just gives an indication of how many
  # elements are in x
  if(is.null(FUN)) {
   cat(".")
   return(x)
  }
  else {
   nfargs<-length(fargs)+1
   farglist<-vector("list",nfargs)
   farglist[[1]]<-x
   if(nfargs > 1) farglist[2:nfargs]<-fargs
   newelement<-do.call(FUN,farglist)
   return(newelement)
  }
 }
}
