traza=function(A)     {
 if (!is.fdata(A)) sum(diag(A),na.rm=TRUE)
 else sum(diag(A$data),na.rm=TRUE)
  }
