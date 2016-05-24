CheckPlotFDRsFormat <- function(lpcfdr.out,frac){
  if(frac<0 || frac>1) stop("frac must be between 0 and 1....")
  if(is.null(lpcfdr.out$fdrlpc) || is.null(lpcfdr.out$fdrt) || length(lpcfdr.out$fdrt)!=length(lpcfdr.out$fdrlpc)) print("lpcfdr.out must have attributes fdrlpc and fdrt, which have equal length....")
} 

