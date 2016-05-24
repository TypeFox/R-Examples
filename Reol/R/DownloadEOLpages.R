DownloadEOLpages <- function(pages, to.file=TRUE, MyKey=NULL, verbose=TRUE, ...) {
  if(Sys.getlocale("LC_ALL") == "C")
    warning("Sys.getlocale is set to C. In order to read UTF characters, you need to set the locale aspect to UTF-8 using Sys.setlocale")
  EOLpages <- vector("list", length=length(pages))
  for (i in sequence(length(pages))) {
    pageNum <- pages[i]
    web <- paste("http://eol.org/api/pages/", pageNum, ".xml?images=75&amp;videos=75&amp;sounds=75&amp;maps=75&amp;text=75&amp;iucn=true&amp;subjects=all&amp;details=true&amp;common_names=true&amp;references=true", sep="")	
	if(!is.null(MyKey))
      web <- paste(web, "&amp;key=", MyKey, sep="")
    if(to.file) {
      write(getURL(web, ...), file=paste("eol", pages[i], ".xml", sep=""))
      if(verbose)
        print(paste("Downloaded ", "eol", pages[i], ".xml", sep=""))
    }
    else {
      EOLpages[[i]] <- getURL(web, ...)
      names(EOLpages)[[i]] <- paste("eol", pages[i], sep="")
      if(verbose)
        print(paste("eol", pages[i], " saved as R object", sep=""))
    }  
    Sys.sleep(1)
  }
  if(to.file)
    return(paste("eol", pages, ".xml", sep=""))
  else
    return(EOLpages)
}
