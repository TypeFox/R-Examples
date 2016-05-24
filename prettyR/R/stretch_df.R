stretch_df<-function(x,idvar,to.stretch,ordervar=NA,include.ordervar=TRUE) {
 xdim <- dim(x)
 # make sure that the fields have names
 if(is.null(names(x))) names(x)<-paste("V",1:xdim[2],sep="")
 # get the names of the fields that are constant within cases
 if(is.na(ordervar))
  nostretch<-names(x)[!(names(x) %in% to.stretch)]
 else nostretch<-names(x)[!(names(x) %in% c(to.stretch,ordervar))]
 nnostretch<-length(nostretch)
  # get the maximum number of records for any case
 maxstretch<-max(table(x[,idvar]))
 # get the indices of the unique IDs
 nsindx<-which(!duplicated(x[,idvar]))
 # get the actual IDs
 IDs<-x[nsindx,idvar]
 # start the new data frame with just the IDs
 newdf<-data.frame(IDs)
 # if there are other constant fields
 if(nnostretch > 1) {
  # add them to the new data frame
  for(newcol in 2:nnostretch) {
   newdf[,newcol]<-x[nsindx,nostretch[newcol]]
  }
 }
 start<-nnostretch
 if(include.ordervar && !is.na(ordervar)) to.stretch<-c(to.stretch,ordervar)
 nstretch<-length(to.stretch)
 for(stretchcol in 1:nstretch) {
  # beginning at the first field beyond the constant fields
  # add a column of NAs 
  for(newcol in (start+1):(start+maxstretch)) newdf[[newcol]]<-NA
  # jump to the next stretched field
  start<-start+maxstretch
 }
 names(newdf)<-c(nostretch,paste(rep(to.stretch,each=maxstretch),
  rep(1:maxstretch,nstretch),sep="_"))
 for(idno in 1:length(IDs)) {
  rows<-which(x[,idvar] == IDs[idno])
  nrows<-length(rows)
  start<-nnostretch
  for(stretchvar in 1:nstretch) {
   if(is.na(ordervar)) stretchorder<-1:nrows
   else stretchorder<-order(x[rows,ordervar])
   newdf[idno,(start+1):(start+nrows)]<-
    x[rows,to.stretch[stretchvar]][stretchorder]
   start<-start+maxstretch
  }
 }
 return(newdf)
}
