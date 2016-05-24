`to.shape`<-function(M,determ=NULL,trace=NULL,first=NULL)
{
 if(all(is.null(determ),is.null(trace),is.null(first))) 
  return(M/det(M)^(1/dim(M)[2]))
 if(is.numeric(determ))
  return(M*(determ/det(M))^(1/dim(M)[2]))
 if(is.numeric(trace))
  return(M*trace/sum(diag(M)))
 M*first/M[1,1]
}

