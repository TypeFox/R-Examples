print.gp=function(x,...){
# default print statement for a gp object
  cat("A Gaussian process object, approximated in the Fourier basis\n")
  cat("Dimension of domain: ",x$d,"\n")
  if(x$d==1){
    cat("Regular grid of size: ",x$gridsize[1],"\n")
  } else{
    cat("Regular grid of size: ",x$gridsize[1],"by ",x$gridsize[2],"\n")
  }
  cat("Spectral density function of correlation function: ",x$specdens.name,"\n")
  cat("Correlation function (spectral density) parameters: ",x$specdens.param,"\n")
}

