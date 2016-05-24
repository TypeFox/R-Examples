rep_n_stack<-function(data,to.stack,stack.names=NULL) {
 datadim<-dim(data)
 stackdim<-dim(to.stack)
 if(is.null(dim(to.stack))) {
  stackdim<-c(1,length(to.stack))
  to.stack<-matrix(to.stack,ncol=length(to.stack))
 }
 # if column names are passed, convert them to indices
 if(is.character(to.stack))
  to.stack<-matrix(which(names(data) %in% to.stack),nrow=stackdim[1],
   ncol=stackdim[2],byrow=TRUE)
 # get the indices of the columns to be replicated
 to.rep<-which(!(1:datadim[2] %in% to.stack))
 nrep<-length(to.rep)
 # if any columns are to be replicated, 
 # form the new data frame by replicating the first to.rep column
 if(nrep) {
  newDF<-data.frame(rep(data[,to.rep[1]],stackdim[2]))
  # if there is more than one to.rep column, add the rest
  if(nrep > 1) {
   for(repvar in 2:nrep)
    newDF[[repvar]]<-rep(data[[to.rep[repvar]]],stackdim[2])
  }
  if(is.null(stackdim)) stackdim<-1
  for(stackrow in 1:stackdim[1]) {
   newDF<-cbind(newDF,rep(names(data[,to.stack[stackrow,]]),each=datadim[1]),
    unlist(data[,to.stack[stackrow,]]))
  }
  if(is.null(stack.names))
   stack.names<-paste(rep(c("group","value"),stackdim[1]),
    rep(1:stackdim[1],each=stackdim[1]),sep="")
  names(newDF)<-c(names(data[to.rep]),stack.names)
 }
 else {
  # if no columns are to be replicated, just string out the whole
  # data frame into one column, then add the column and row names
  if(is.data.frame(data)) {
   newDF<-data.frame(unlist(data))
   newDF[[2]]<-rep(names(data),datadim[2])
  }
  else {
   newDF<-data.frame(as.vector(data))
   newDF[[2]]<-rep(colnames(data),datadim[2])
  }
  newDF[[3]]<-rep(rownames(data),each=datadim[1])
  if(!is.null(stack.names)) names(newDF)<-stack.names
 }
 rownames(newDF)<-NULL
 return(newDF)
}
