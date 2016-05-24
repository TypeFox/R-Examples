easyFeatures <-
function()
{
	feaStatWindowing<-.easyFeaStatWindowing()
	feaDoubleStatWindowing<-.easyFeaDoubleStatWindowing()
	feaSpecStatWindowing<-.easyFeaSpecStatWindowing()
	feaPCA<-.easyFeaPCA()
	feaCWT<-.easyFeaCWT()

	result<-list(feaStatWindowing=feaStatWindowing,feaDoubleStatWindowing=feaDoubleStatWindowing,
	feaSpecStatWindowing=feaSpecStatWindowing, feaPCA=feaPCA,
	feaCWT=feaCWT)

	class(result)<-"Features"
	return(result)
}


print.Features <-
  function(x,...){
    cat("Features:")
    cat("\n\nWindowing of the data. \nNumber of applications: ",x$feaStatWindowing$N)
    cat("\n\nDouble windowing of the data. \nNumber of applications: ",x$feaDoubleStatWindowing$N)
    cat("\n\nWindowing of the spectrum. \nNumber of applications: ",x$feaSpecStatWindowing$N)
    if(x$feaPCA$feaLoadingsPCA)
    {
      cat("\n\nLoadings of PCA: Yes")
    }else{
      cat("\n\nLoadings of PCA: No")
    }
    
    if(x$feaPCA$feaSpecLoadingsPCA)
    {
      cat("\n\nLoadings of PCA over the spectrum: Yes")
    }else{
      cat("\n\nLoadings of PCA over the spectrum: No")
    }
    cat("\n\nSignals of principal components of PCA. \nNumber of applications: ",x$feaPCA$feaSignalsPCA)
    
    
    if(x$feaCWT$feaCWT)
    {
      cat("\n\nContinuous Wavelet Transform: Yes\n")
    }else{
      cat("\n\nContinuous Wavelet Transform: No\n")
    }
  }


summary.Features <-
  function(object, ...){
    x<- object
    cat("Features:")
    cat("\n\nWindowing of the data. \nNumber of applications: ",x$feaStatWindowing$N)
    cat("\n\nDouble windowing of the data. \nNumber of applications: ",x$feaDoubleStatWindowing$N)
    cat("\n\nWindowing of the spectrum. \nNumber of applications: ",x$feaSpecStatWindowing$N)
    if(x$feaPCA$feaLoadingsPCA)
    {
      cat("\n\nLoadings of PCA: Yes")
    }else{
      cat("\n\nLoadings of PCA: No")
    }
    
    if(x$feaPCA$feaSpecLoadingsPCA)
    {
      cat("\n\nLoadings of PCA over the spectrum: Yes")
    }else{
      cat("\n\nLoadings of PCA over the spectrum: No")
    }
    cat("\n\nSignals of principal components of PCA. \nNumber of applications: ",x$feaPCA$feaSignalsPCA)
    
    
    if(x$feaCWT$feaCWT)
    {
      cat("\n\nContinuous Wavelet Transform: Yes\n")
    }else{
      cat("\n\nContinuous Wavelet Transform: No\n")
    }
  }
