# Arguments:
# x - an R object. If a list, it will be searched for non-list elements
# FUN - an optional function to perform on the non-list elements of x
# fargs - additional arguments to FUN (removed until I can get it to work)
# max - whether to search for minimal or maximal values of FUN (or elements of x)
# retval - a list of three components (see below)
# Return value:
# The final value of the list retval
#  indx - a vector of the list indices of the element of x with the extreme value
#  element - the element of x with the extreme value
#  value - the extreme value of the first element of x or of FUN applied to x
# listCrawler descends an R list until it encounters a non-list element.
# it then returns either the first value of that element, or if FUN is not NULL,
# the result of applying FUN to the element plus the element itself

listCrawler<-function(x,FUN=NULL,maxval=TRUE,
 retval=list(indx=vector("numeric",0),element=NULL,value=NA)) {

 if(is.list(x)) {
  thisindx<-retindx<-1
  for(lindex in 1:length(x)) {
   # recursively run through the elements of this list
   newretval<-listCrawler(x[[lindex]],FUN=FUN,maxval=maxval,retval=retval)
   # grab the first element and value in case the latter is the extreme
   if(is.na(retval$value[1])) {
    retval$element<-newretval$element
    retval$value<-newretval$value
   }
   # maxval signals whether to look for a max or min in the first element of
   # the value returned as it might be a vector
   selectcond<-
    ifelse(maxval,newretval$value[1]>retval$value[1],newretval$value[1]<retval$value[1])
   if(selectcond) {
    # insert the new extreme value and element and 
    # grab the current index and lower ones
    retval$value<-newretval$value
    retval$element<-newretval$element
    retindx<-newretval$indx
    thisindx<-lindex
   }
  }
  # now concatenate the current index for the extremum and lower ones
  retval$indx<-c(thisindx,retindx)
  return(retval)
 }
 else {
  if(is.null(FUN)) {
   # just use the first element of x as value and x as element
   retval$value<-x[1]
   retval$element<-x
  }
  else {
   # if there is something to do, do it
   retval$value<-do.call(FUN,list(x))
   retval$element<-x
  }
 }
 return(retval)
}
