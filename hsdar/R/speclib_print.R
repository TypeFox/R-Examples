setMethod("print", signature(x = "Speclib"), 
          function(x)
{
  cat(paste("Summary of ", class(x), "\n", sep = ""))
  .printUsagehistory(x)
  cat("\n\nSummary of spectra\n")
  cat("---------------------\n")
  cat(paste("Total number of spectra :",dim(x)[1]))
  cat(paste("\nNumber of bands :",dim(x)[2]))
  if (length(x@fwhm)==1)
  {
    cat(paste("\nWidth of bands :",x@fwhm))
  } else {
    cat(paste("\nMean width of bands :",round(mean(x@fwhm),2)))
  }
  cat("\n")
  if (x@spectra@fromRaster)
  {
    cat(paste("Use RasterBrick for spectra"))
    if (x@spectra@spectra_ra@data@inmemory)
    {
      cat(paste(" (in RAM)\n"))
    } else {
      cat(paste(" stored at\n'", x@spectra@spectra_ra@file@name, "'\n", sep =""))
    }
  }    
}
)

setMethod ('show' , signature(object = "Speclib"), 
           function(object)
{
  print(object)  
}
)

.printUsagehistory <- function(x)
{
  if (length(x@usagehistory)>0)
  {
    cat("\n\nHistory of usage\n")
    cat("---------------------\n")
    for (i in 1:length(x@usagehistory)) cat(paste("(",i,")   ",x@usagehistory[i],"\n",sep=""))
  }
}