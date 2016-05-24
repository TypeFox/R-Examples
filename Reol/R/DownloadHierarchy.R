ProviderCount <- function(MyEOLs, verbose=FALSE) {
  #Returns vector of provider coverage
  results <- GatherProviderDataFrame(MyEOLs, extended.output=FALSE)[,-1:-2]
  results <- results[, -dim(results)[2]]
  counts <- apply(results, 2, sum)
  names(counts) <- colnames(results)
  counts <- sort(counts, decreasing=TRUE)
  if (verbose)
    print(t(t(counts)))	
  return(counts)
}

BestProvider <- function(MyEOLs) {
  #Returns the provider with the most taxonomic coverage
  return(names(ProviderCount(MyEOLs))[1])	
}


DownloadHierarchy <- function(MyEOLs, to.file=TRUE, database=NULL, verbose=TRUE, ...) {
#MyEOLs can be a file or an R object
#to.file is whether you want to save the information as a file (T) or an R object (F)
  #Downloads provider database
  if(is.null(database))
    database <- BestProvider(MyEOLs)	
  results <- GatherProviderDataFrame(MyEOLs, extended.output=TRUE)
  column <- which(colnames(results) == paste(database, ".taxonID", sep=""))
  pages <- results[,column] 
  hierpages <- vector("list", length=length(pages))
  for (i in sequence(length(pages))) {
    if (!is.na(pages[i])) {
      pageNum<-pages[i]
      web <- paste("http://eol.org/api/hierarchy_entries/1.0/", pageNum, sep="")
      if(to.file) {
        write(getURL(web, ...), file=paste("hier", pages[i], ".xml", sep=""))
        if(verbose)
          print(paste("Downloaded ", "hier", pages[i], ".xml", sep=""))
      }
      else {
        hierpages[[i]] <- getURL(web, ...)
        names(hierpages)[[i]] <- paste("hier", pages[i], sep="")
        if(verbose)
          print(paste("hier", pages[i], " saved as R object", sep=""))
      }
      Sys.sleep(1)
    }
  }
  if(to.file)
    return(paste("hier", pages, ".xml", sep=""))
  else
    return(hierpages)
}
