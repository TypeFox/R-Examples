`makePathName` <- function# create appropriate directory structure if needed
### create appropriate directory structure if needed
(path, ##<< path to create
 MakePath=TRUE, ##<< if yes, create directory if not exists
 verbose=0 ##<< level of verbosity
 ){
  N <- nchar(path);
  if (substring(path,N,N) != "/"){
  	path <- paste(path, "/",sep="")
  }
   if (MakePath)
   if (!file.exists(path)) {
  	tmp <- dir.create(path);
  	if (!tmp) cat(path, "could not be created.\n");
  	if (tmp & verbose) cat(path, "successfully created.\n");
  	flush(stdout())
  }

  return(path)
### returns absolute path
}

