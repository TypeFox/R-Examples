Nrows <- function(ndm,nvrr,ncon2,rtheta2){
  # res holds the results for the fortran call 
  res <- .Fortran("cnvrt", as.integer(ndm),as.integer(nvrr),as.integer(ncon2),as.single(rtheta2),
                  nvrtr=integer(1),rxvt=single(45000),ifa=integer(4))
  # returns the number of rows in the output matrix rxvt
     return(res$nvrtr)
  }