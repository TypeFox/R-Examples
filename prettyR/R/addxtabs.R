addxtabs<-function(x) {
 if(!is.list(x) || class(x[[1]]) != "xtab")
  stop("addxtabs adds a list of xtab results with equal numbers of columns")
 allrownames<-unlist(lapply(lapply(x,"[[","counts"),rownames))
 rowlabels<-sort(unique(allrownames))
 allcolnames<-unlist(lapply(lapply(x,"[[","counts"),colnames))
 collabels<-sort(unique(allcolnames))
 dimx1<-dim(x[[1]]$counts)
 outmat<-matrix(0,nrow=length(rowlabels),ncol=dimx1[2])
 for(xt in 1:length(x)) {
  currowlab<-rownames(x[[xt]]$counts)
  for(row in 1:length(currowlab)) {
   curcollab<-colnames(x[[xt]]$counts)
   whichrow<-rowlabels%in%currowlab[row]
   for(column in 1:length(curcollab)) {
    whichcol<-collabels%in%curcollab[column]
    outmat[whichrow,whichcol]<-outmat[whichrow,whichcol]+x[[xt]]$counts[row,column]
   }
  }
 }
 rownames(outmat)<-rowlabels
 colnames(outmat)<-collabels
 return(outmat)
}
