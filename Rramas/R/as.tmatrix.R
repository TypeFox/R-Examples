as.tmatrix <-
function(x, names.st=NULL,...){
   x <- as.matrix(x,...)
   if(diff(dim(x)) !=0) stop("only square matrices can be considered
                              as a transition matrix.")
   di <- dim(x)[[1]]
   m.names <- dimnames(x)[[1]]
   if(is.null(m.names)) m.names <- names.st
   if(is.null(m.names)) m.names <- paste("stage.",1:di ,sep="")
   
   dimnames(x) <- list(m.names, m.names)
   class(x)<-c("tmatrix", class(x))
  return(x)
 }

