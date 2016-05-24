# muranks substitutes the mean of unallocated ranks for missing
# values in a set of ranks before calculating the mean rankings for data objects.
# This procedure assumes that the ranking method (usually a respondent or judge)
# assigns ranks somewhere within a specified range 1:n.
# Any ranks within this range may be assigned.
# The data are assumed to look like:
# NA  1   3   2   NA  4   NA  9.5  9.5  NA
# which is a possible ranking of ten data objects where six ranks with two
# ties have been allocated. Note that all missing ranks (5,6,7,8) are adjacent.
# If this is not the case in a row, the function will delete that row from
# the vector or matrix that is returned. This prevents counterintuitive rank
# imputations like 1 2 4 5 NA NA -> 1 2 4 5 4.5 4.5
# However, for "competition ranks" see below.
# Tied ranks will be expanded to exclude all ranks that would
# have been present if the ties had not occurred in the calculation of the
# mean of unallocated ranks. Thus in the example above, if expand.ties is TRUE,
# only ranks 3,6,9 and 10 are considered missing, and unallocated ranks are
# replaced with 7. This method is probably what users will want if ties are
# given the mean of the tied rankings as in the default method of the rank
# function. If ties are not scored as the mean rank, this method will _not_ 
# work properly.
# To process "competition ranks", set rankx to TRUE.
# Be careful, this will process _any_ set of numbers.

muranks<-function(x,allranks=NULL,rankx=FALSE) {
 dimx<-dim(x)
 if(is.null(dimx)) {
  if(rankx) x<-rank(x,na.last="keep")
  lenx<-length(x)
  # assume that ranks are from 1 to the number of rankings
  if(is.null(allranks)) allranks<-1:lenx
  # discard any vector with a value outside the range of ranks
  if(any(!x[!is.na(x)]%in%allranks)) return(NULL)
  mr<-allranks[!allranks%in%x]
  # discard any vector with non-sequential missing ranks
  if(any(diff(mr[order(mr)])>1) | 1%in%mr) return(NULL)
  xx<-x
  # perform the expansion on a duplicate matrix to maintain the ties.
  for(i in 1:lenx) {
   nties<-sum(xx==xx[i],na.rm=TRUE)
   if(nties > 1) {
    tiesdiv<-ifelse(xx[i]==floor(xx[i]),3,4)
    xx[which(xx==xx[i])]<-(xx[i]-nties/tiesdiv):(xx[i]+nties/tiesdiv)
   }
  }
  subrank<-mean(allranks[!allranks%in%xx])
  x[which(is.na(x))]<-subrank
 }
 else {
  if(is.null(allranks))
   allranks<-1:max(c(max(x,na.rm=TRUE),dimx[2]),na.rm=TRUE)
  if(is.data.frame(x)) x<-as.matrix(x)
  badrows<-rep(FALSE,dimx[1])
  for(i in 1:dimx[1]) {
   newrow<-muranks(x[i,],allranks=allranks,rankx=rankx)
   if(is.null(newrow)) badrows[i]<-TRUE
   else x[i,]<-newrow
  }
  if(any(badrows)) {
   cat(sum(badrows),"rows were discarded due to invalid rank values\n")
   x<-x[-which(badrows),]
  }
 }
 return(x)
}
